// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2023 INRIA
  
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::DamarisWriter
 */
#ifndef EWOMS_DAMARIS_WRITER_HH
#define EWOMS_DAMARIS_WRITER_HH

#include <dune/grid/common/partitionset.hh>

#include <ebos/collecttoiorank.hh>
#include <ebos/eclbasevanguard.hh>
#include <ebos/eclgenericwriter.hh>
#include <ebos/ecloutputblackoilmodule.hh>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <limits>
#include <stdexcept>
#include <string>

#include <fmt/format.h>

#include <opm/simulators/utils/GridDataOutput.hpp>
#include <damaris/util/DamarisGeometryData.hpp>


namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutput {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableDamarisOutputCollective {
    using type = UndefinedProperty;
};
} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass them to Damaris server processes.
 *
 * Currently only passing through PRESSURE, GLOBAL_CELL_INDEX and MPI_RANK information.
 * This clss will be enhanced to pass through the 3D mesh information to Damaris to enable
 * in situ visualization via Paraview or Ascent. And developed so that variables specified 
 * through the Eclipse input deck will be available to Damaris.
 */
 
 
template <class TypeTag>
class DamarisWriter : public EclGenericWriter<GetPropType<TypeTag, Properties::Grid>,
                                          GetPropType<TypeTag, Properties::EquilGrid>,
                                          GetPropType<TypeTag, Properties::GridView>,
                                          GetPropType<TypeTag, Properties::ElementMapper>,
                                          GetPropType<TypeTag, Properties::Scalar>>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    
    using BaseType = EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>;
    
public:
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableDamarisOutputCollective,
                             "Write output via Damaris using parallel HDF5 to get single file per timestep instead of one per Damaris core.");
    }

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    DamarisWriter(Simulator& simulator)
        : BaseType(simulator.vanguard().schedule(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().summaryConfig(),
                   simulator.vanguard().grid(),
                   simulator.vanguard().grid().comm().rank() == 0 ? &simulator.vanguard().equilGrid() : nullptr,
                   simulator.vanguard().gridView(),
                   simulator.vanguard().cartesianIndexMapper(),
                   simulator.vanguard().grid().comm().rank() == 0 ? &simulator.vanguard().equilCartesianIndexMapper() : nullptr,
                   false, false)
        , simulator_(simulator)
    {
        this->damarisUpdate_ = true ;
        
        rank_ = simulator_.vanguard().grid().comm().rank() ;
        nranks_ = simulator_.vanguard().grid().comm().size();
        
        const auto& gridView = simulator_.gridView();
        const auto& interior_elements = elements(gridView, Dune::Partitions::interior);
        // Get the size of the unique vector elements (excludes the shared 'ghost' elements)
        numElements_ = std::distance(interior_elements.begin(), interior_elements.end());
        
        this->elements_rank_offsets_.resize(nranks_) ;
        this->damarisOutputModule_ = std::make_unique<EclOutputBlackOilModule<TypeTag>>(simulator, this->collectToIORank_);
    }

    ~DamarisWriter()
    { }

    /*!
     * \brief Writes localCellData through to Damaris servers. Sets up the unstructured mesh which is passed to Damaris.
     */
    void writeOutput(data::Solution& localCellData , bool isSubStep)
    {
        OPM_TIMEBLOCK(writeOutput);
        const int reportStepNum = simulator_.episodeIndex() + 1;
        
        // added this as localCellData was not being written
        if (!isSubStep)
            this->damarisOutputModule_->invalidateLocalData() ;  
        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->damarisOutputModule_->outputErrorLog(simulator_.gridView().comm());

        // The damarisWriter is not outputing well or aquifer data (yet)
        auto localWellData = simulator_.problem().wellModel().wellData(); // data::Well
        
        if (! isSubStep) 
        {
            if (localCellData.size() == 0) {
                this->damarisOutputModule_->assignToSolution(localCellData);
            }

            // add cell data to perforations for Rft output
            this->damarisOutputModule_->addRftDataToWells(localWellData, reportStepNum);
            
            // On first call and if the mesh and variable size change then set damarisUpdate_ to true
            if (damarisUpdate_ == true) {
                // Sets the damaris parameter values "n_elements_local" and "n_elements_total" 
                // which define sizes of the Damaris variables, per-rank and globally (over all ranks).
                // Also sets the offsets to where a ranks array data sits within the global array. 
                // This is usefull for HDF5 output and for defining distributed arrays in Dask.
                this->setupDamarisWritingPars(simulator_.vanguard().grid().comm(), numElements_, elements_rank_offsets_);
                
                // sets data for non-time-varying variables MPI_RANK and GLOBAL_CELL_INDEX
                this->SetGlobalIndexForDamaris() ; 
                
                // Adds all mesh variables (x,y,z, type, offsets and coords) and sets damarisUpdate_  to false
                this->writeDamarisGridOutput(isSubStep) ;
                
                // Currently by default we assume static grid (unchanging through the simulation)
                // Set damarisUpdate_ to true if we want to update the geometry to sent to Damaris 
                this->damarisUpdate_ = false; 
            }
            
            
            if (this->damarisOutputModule_->getPRESSURE_ptr() != nullptr) 
            {
                int64_t temp_int64_t[1];
                temp_int64_t[0] = static_cast<int64_t>(this->elements_rank_offsets_[rank_]);
                dam_err_ = damaris_set_position("PRESSURE", temp_int64_t);
                if (dam_err_ != DAMARIS_OK && rank_ == 0) {
                    OpmLog::error(fmt::format("ERORR: damariswriter::writeOutput()       : ( rank:{}) damaris_set_position(PRESSURE, ...), Damaris Error: {}  ",  rank_, damaris_error_string(dam_err_) ));
                }
                
                dam_err_ = damaris_write("PRESSURE", (void*)this->damarisOutputModule_->getPRESSURE_ptr());
                if (dam_err_ != DAMARIS_OK) {
                   OpmLog::error(fmt::format("ERORR: damariswriter::writeOutput()       : ( rank:{}) damaris_write(PRESSURE, ...), Damaris Error: {}  ",  rank_, damaris_error_string(dam_err_) ));
                }
            }
            
            dam_err_ =  damaris_end_iteration();
            if (dam_err_ != DAMARIS_OK) {
                std::cerr << "ERROR rank =" << rank_ << " : damariswriter::writeOutput() : damaris_end_iteration()" 
                << ", Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
            }
         } // end of ! isSubstep
    }
 
private:

    int dam_err_ ;
    int rank_  ;       
    int nranks_ ;
    int numElements_ ;  ///<  size of the unique vector elements
    
    Simulator& simulator_;
    std::unique_ptr<EclOutputBlackOilModule<TypeTag>> damarisOutputModule_;
    std::vector<unsigned long long> elements_rank_offsets_ ;
    bool damarisUpdate_ = false;  ///< Whenever this is true writeOutput() will set up Damaris mesh information and offsets of model fields

    static bool enableDamarisOutput_()
    { 
        return EWOMS_GET_PARAM(TypeTag, bool, EnableDamarisOutput); 
    }

    void SetGlobalIndexForDamaris () 
    {
        // GLOBAL_CELL_INDEX is used to reorder variable data when writing to disk 
        // This is enabled using select-file="GLOBAL_CELL_INDEX" in the <variable> XML tag
        if ( this->collectToIORank_.isParallel() ){
            const std::vector<int>& local_to_global =  this->collectToIORank_.localIdxToGlobalIdxMapping(); 
            dam_err_ = damaris_write("GLOBAL_CELL_INDEX", local_to_global.data());
        } else {
            std::vector<int> local_to_global_filled ;
            local_to_global_filled.resize(this->numElements_) ;
            for (int i = 0 ; i < this->numElements_ ; i++)
            {
                local_to_global_filled[i] = i ;
            }
            dam_err_ = damaris_write("GLOBAL_CELL_INDEX", local_to_global_filled.data());
        }

        if (dam_err_ != DAMARIS_OK) {
            std::cerr << "ERROR rank =" << rank_ << " : eclwrite::writeOutput() : damaris_write(\"GLOBAL_CELL_INDEX\", local_to_global.data())" 
                  << ", Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
        }

        std::vector<int> mpiRank(this->numElements_, rank_ ) ;
        dam_err_ = damaris_write("MPI_RANK", mpiRank.data() ) ;
        if (dam_err_ != DAMARIS_OK) {
           std::cerr << "ERROR rank =" << rank_ << " : damariswriter::writeOutput() : damaris_write(\"MPI_RANK\""
           << "...), failed. Damaris error = " <<  damaris_error_string(dam_err_) << std::endl ;
        }
           
    }

    void setupDamarisWritingPars(Parallel::Communication comm, const int n_elements_local_grid, std::vector<unsigned long long>& elements_rank_offsets)
    {
        const int nranks = comm.size();
        const int rank = comm.rank();

        std::vector<unsigned long long> elements_rank_sizes(nranks); // one for each rank -- to be gathered from each client rank
        // n_elements_local_grid should be the full model size
        const unsigned long long n_elements_local = n_elements_local_grid;

        // This gets the n_elements_local from all ranks and copies them to a std::vector of all the values on all ranks
        // (elements_rank_sizes[]).
        comm.allgather(&n_elements_local, 1, elements_rank_sizes.data());
        elements_rank_offsets[0] = 0ULL;
        // This scan makes the offsets to the start of each ranks grid section if each local grid data was concatenated (in
        // rank order)
        for (int t1 = 1; t1 < nranks; t1++) {
            elements_rank_offsets[t1] = elements_rank_offsets[t1 - 1] + elements_rank_sizes[t1 - 1];
        }

        // find the global/total size
        unsigned long long n_elements_global_max = elements_rank_offsets[nranks - 1];
        n_elements_global_max += elements_rank_sizes[nranks - 1]; // add the last ranks size to the already accumulated offset values

        if (rank == 0) {
            OpmLog::debug(fmt::format("In setupDamarisWritingPars(): n_elements_global_max = {}", n_elements_global_max));
        }

        // Set the paramater so that the Damaris servers can allocate the correct amount of memory for the variabe
        // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
        // ToDo: Do we need to check that local ranks are 0 based ?
        int temp_int = static_cast<int>(elements_rank_sizes[rank]);
        dam_err_ = damaris_parameter_set("n_elements_local", &temp_int, sizeof(int));
        if (dam_err_ != DAMARIS_OK && rank == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_parameter_set(\"n_elements_local\", &temp_int, sizeof(int));");
        }
        // Damaris parameters only support int data types. This will limit models to be under size of 2^32-1 elements
        // ToDo: Do we need to check that n_elements_global_max will fit in a C int type (INT_MAX)
        temp_int = static_cast<int>(n_elements_global_max);
        dam_err_ = damaris_parameter_set("n_elements_total", &temp_int, sizeof(int));
        if (dam_err_ != DAMARIS_OK && rank == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_parameter_set(\"n_elements_total\", &temp_int, sizeof(int));");
        }

        // Use damaris_set_position to set the offset in the global size of the array.
        // This is used so that output functionality (e.g. HDF5Store) knows global offsets of the data of the ranks
        int64_t temp_int64_t[1];
        temp_int64_t[0] = static_cast<int64_t>(elements_rank_offsets[rank]);
        dam_err_ = damaris_set_position("PRESSURE", temp_int64_t);
        if (dam_err_ != DAMARIS_OK && rank == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_set_position(\"PRESSURE\", temp_int64_t);");
        }
        dam_err_ = damaris_set_position("GLOBAL_CELL_INDEX", temp_int64_t);
        if (dam_err_ != DAMARIS_OK && rank == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_set_position(\"GLOBAL_CELL_INDEX\", temp_int64_t);");
        }
        dam_err_ = damaris_set_position("MPI_RANK", temp_int64_t);
        if (dam_err_ != DAMARIS_OK && rank == 0) {
            OpmLog::error("Damaris library produced an error result for "
                          "damaris_set_position(\"PRESSURE\", temp_int64_t);");
        }
    }


    void writeDamarisGridOutput(bool isSubStep)
    {
        const auto& gridView = simulator_.gridView();
        //  const auto& interior_elements = elements(gridView, Dune::Partitions::interior);
        //  Get the size of the unique vector elements (excludes the shared 'ghost' elements)
        //  const int numElements = std::distance(interior_elements.begin(), interior_elements.end());
        //  Sets the damaris parameter values which then defines sizes of the arrays per-rank and globally.
        //  Also sets the offsets to where a ranks array data sits within the global array. 
        //  This is usefull for HDF5 output and for defining distributed arrays in Dask.
        //  Opm::DamarisOutput::setupDamarisWritingPars(simulator_.vanguard().grid().comm(), numElements);
        
        Opm::GridDataOutput::SimMeshDataAccessor geomData(gridView, Dune::Partitions::interior) ; // N.B. we cannot reuse the same object using a different partition.
        try {
           
            damaris::model::vertex_data_structure vertex_structure = damaris::model::VERTEX_SEPARATE_X_Y_Z ;  // define this as we know it works with Ascent
            damaris::model::DamarisGeometryData damarisMeshVars(vertex_structure, geomData.getNVertices(), 
                                                        geomData.getNCells(), geomData.getNCorners(), rank_) ;
           
            const bool hasPolyCells = geomData.polyhedralCellPresent() ;
            if ( hasPolyCells ) {
                std::cout << "The DUNE geometry grid has polyhedral elements - currently not supported by Damaris " << std::endl ;
            } 
            damarisMeshVars.set_hasPolyhedralCells(hasPolyCells) ;
           
            // This is our template XML model for x,y,z coordinates
            // <parameter name="n_coords_local"     type="int" value="1" />
            // <parameter name="n_coords_global"    type="int" value="1" comment="only needed if we need to write to HDF5 in Collective mode/>
            // <layout    name="n_coords_layout"    type="double" dimensions="n_coords_local"   comment="For the individual x, y and z coordinates of the mesh vertices, these values are referenced in the topologies/topo/subelements/connectivity_pg data"  />
            // <group name="coordset/coords/values"> 
            //     <variable name="x"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
            //     <variable name="y"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
            //     <variable name="z"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonConduitTest" time-varying="false" />
            // </group>
            // 
            int xyz_coord_dims = 1 ;
            std::vector<std::string> param_names = {"n_coords_local"} ;  // a vector of strings as a variables layout may be defined by multiple parameters 
            std::string variable_x = "coordset/coords/values/x" ;  // This string must match the group/variable name of the Damaris XML file 
            std::string variable_y = "coordset/coords/values/y" ;  // This string must match the group/variable name of the Damaris XML file 
            std::string variable_z = "coordset/coords/values/z" ;  // This string must match the group/variable name of the Damaris XML file 
            
            damarisMeshVars.set_damaris_var_name_vertex_x(xyz_coord_dims,  param_names, variable_x) ;
            damarisMeshVars.set_damaris_var_name_vertex_y(xyz_coord_dims,  param_names, variable_y) ;
            damarisMeshVars.set_damaris_var_name_vertex_z(xyz_coord_dims,  param_names, variable_z) ;
            
            // Used to store Damaris parameters, which are then used to resize the 
            // shared memory region used to save the mesh data
            // there should be as many values for paramaters as names used in param_names.
            // Each name corresponds to the value at the same position in the vector.
            std::vector<int> param_vertices ;  
            std::vector<int> param_connectivity ;
            std::vector<int> param_offsets ;
            
            
            param_vertices.push_back(geomData.getNVertices() ) ; 
            // For the vertex data x, y and z arrays, SetAll_VERTEX_SEPARATE_X_Y_Z_shmem()  will set the arrays 
            // size (using the paramater value) and then allocate the shared memory region. This is where we 
            // will write the vertex data to, so that Damaris has access to it.
            damarisMeshVars.SetAll_VERTEX_SEPARATE_X_Y_Z_shmem(param_vertices) ;  
            
            // Now we can return the memory that Damaris has allocated in shmem 
            damaris::model::DamarisVar<double>* var_x =  dynamic_cast<damaris::model::DamarisVar<double>* >(damarisMeshVars.get_x()) ;
            damaris::model::DamarisVar<double>* var_y =  dynamic_cast<damaris::model::DamarisVar<double>* >(damarisMeshVars.get_y()) ;
            damaris::model::DamarisVar<double>* var_z =  dynamic_cast<damaris::model::DamarisVar<double>* >(damarisMeshVars.get_z()) ;
            
            if ( geomData.writeGridPoints(var_x->data_ptr(),var_y->data_ptr(),var_z->data_ptr()) < 0)
                 DUNE_THROW(Dune::IOError, geomData.getError()  );
            
           
            // We do not need these as the ~DamarisGeometryData destructor will call them  
            // damarisMeshVars.CommitAll_VERTEX_SEPARATE_X_Y_Z_shmem() ;
            // damarisMeshVars.ClearAll_VERTEX_SEPARATE_X_Y_Z_shmem() ;
            
            //  This is our template XML model for connectivity 
            // <parameter name="n_connectivity_ph"        type="int"  value="1" />
            // <layout    name="n_connections_layout_ph"  type="int"  dimensions="n_connectivity_ph"   comment="Layout for connectivities "  />
            // <parameter name="n_offsets_types_ph"       type="int"  value="1" />
            // <layout    name="n_offsets_layout_ph"      type="int"  dimensions="n_offsets_types_ph"  comment="Layout for the offsets_ph"  />
            // <layout    name="n_types_layout_ph"        type="char" dimensions="n_offsets_types_ph"  comment="Layout for the types_ph "  />
            // <group name="topologies/topo/elements">
            //     <variable name="connectivity" layout="n_connections_layout_ph"  type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
            //     <variable name="offsets"      layout="n_offsets_layout_ph"    type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
            //     <variable name="types"        layout="n_types_layout_ph"    type="scalar"  visualizable="false"  unit=""   script="PythonConduitTest" time-varying="false" />
            // </group>
            // 
           
            param_names[0] = "n_connectivity_ph" ; 
            std::string varname = std::string("topologies/topo/elements/connectivity") ;  // This string must match the group/variable name of the Damaris XML file 
            damarisMeshVars.set_damaris_var_name_connectivity(xyz_coord_dims, param_names , varname ) ;
            param_names[0] = "n_offsets_types_ph" ; 
            varname = std::string("topologies/topo/elements/offsets") ;                   // This string must match the group/variable name of the Damaris XML file 
            damarisMeshVars.set_damaris_var_name_offsets(xyz_coord_dims, param_names, varname ) ;
            param_names[0] = "n_offsets_types_ph" ; 
            varname = std::string("topologies/topo/elements/types") ;                     // This string must match the group/variable name of the Damaris XML file 
            damarisMeshVars.set_damaris_var_name_types(xyz_coord_dims, param_names, varname ) ;
            
            // Here we retrieve the DamarisVar objects. We need to *match the type* of the data as set in the Damaris XML file
            damaris::model::DamarisVar<int>* var_connectivity =  dynamic_cast<damaris::model::DamarisVar<int>* >(damarisMeshVars.get_connectivity()) ;
            damaris::model::DamarisVar<int>* var_offsets      =  dynamic_cast<damaris::model::DamarisVar<int>* >(damarisMeshVars.get_offsets()) ;
            damaris::model::DamarisVar<char>* var_types       =  dynamic_cast<damaris::model::DamarisVar<char>* >(damarisMeshVars.get_types()) ;
            
            // Set the Damaris shared memory for each variable needed to store mesh data
            param_connectivity.push_back( geomData.getNCorners() ) ;
            var_connectivity->SetDamarisParameter(param_connectivity) ;
            var_connectivity->SetPointersToDamarisShmem() ;
           
            param_offsets.push_back(geomData.getNCells() ) ;
            var_offsets->SetDamarisParameter(param_offsets) ;
            var_offsets->SetPointersToDamarisShmem() ;
            var_types->SetDamarisParameter(param_offsets) ;
            var_types->SetPointersToDamarisShmem() ;
            
            // Copy the mesh data from the Durne grid
            long i = 0 ;
            Opm::GridDataOutput::ConnectivityVertexOrder vtkorder = Opm::GridDataOutput::VTK ;
            
            i = geomData.writeConnectivity(var_connectivity->data_ptr(), vtkorder) ;
            if ( i  != geomData.getNCorners())
                 DUNE_THROW(Dune::IOError, geomData.getError() );
            
            i = geomData.writeOffsetsCells(var_offsets->data_ptr()) ;
            if ( i != geomData.getNCells()+1)
                 DUNE_THROW(Dune::IOError,geomData.getError() );
            
            i = geomData.writeCellTypes(var_types->data_ptr()) ;
            if ( i != geomData.getNCells())
                 DUNE_THROW(Dune::IOError,geomData.getError() );
            //  Commit and clear damaris functions are called when the object goes out of scope (in the destructor)
        }
        catch (std::exception& e) 
        {
            std :: cout << e.what() << std::endl;
        }
    }

    void prepareLocalCellData(const bool isSubStep,
                              const int  reportStepNum)
    {
        OPM_TIMEBLOCK(prepareLocalCellData);
        if (damarisOutputModule_->localDataValid()) {
            return;
        }

        const auto& gridView = simulator_.vanguard().gridView();
        const int numElements = gridView.size(/*codim=*/0);
        const bool log = this->collectToIORank_.isIORank();

        damarisOutputModule_->allocBuffers(numElements, reportStepNum,
                                      isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        {
        OPM_TIMEBLOCK(prepareCellBasedData);
        for (const auto& elem : elements(gridView)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            damarisOutputModule_->processElement(elemCtx);
        }
        }
        if(!simulator_.model().linearizer().getFlowsInfo().empty()){
            OPM_TIMEBLOCK(prepareFlowsData);
            for (const auto& elem : elements(gridView)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                damarisOutputModule_->processElementFlows(elemCtx);
            }
        }
        {
        OPM_TIMEBLOCK(prepareBlockData);
        for (const auto& elem : elements(gridView)) {
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            damarisOutputModule_->processElementBlockData(elemCtx);
        }
        }
        {
        OPM_TIMEBLOCK(prepareFluidInPlace);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int dofIdx=0; dofIdx < numElements; ++dofIdx){
                const auto& intQuants = *(simulator_.model().cachedIntensiveQuantities(dofIdx, /*timeIdx=*/0));
                const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);
                damarisOutputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
        }
        }
        damarisOutputModule_->validateLocalData();
        OPM_END_PARALLEL_TRY_CATCH("DamarisWriter::prepareLocalCellData() failed: ", simulator_.vanguard().grid().comm());
    }

};
} // namespace Opm

#endif
