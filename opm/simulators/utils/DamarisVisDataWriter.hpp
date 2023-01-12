// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef OPM_DAMARIS_DATA_WRITER_HH
#define OPM_DAMARIS_DATA_WRITER_HH

#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <type_traits>
#include <vector>
#include <list>
#include <map>

#include <dune/common/visibility.hh>
#include <dune/common/typetraits.hh>
// #include <dune/common/exceptions.hh>
#include <dune/common/indent.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/path.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/dataarraywriter.hh>
//#include <dune/grid/io/file/vtk/function.hh>
//#include <dune/grid/io/file/vtk/pvtuwriter.hh>
//#include <dune/grid/io/file/vtk/streams.hh>
//#include <dune/grid/io/file/vtk/vtuwriter.hh>

#include <Damaris.h>

/** @file
    @author Joshua Bowden
    @brief Alows model geometry specification to be passed to Damaris. Based off the Dune VTKWriter.hh code
 */

/**
  // From the opm-simulators repository
  #include <opm/simulators/utils/DamarisVisDataWriter.hpp>
  ...
  Dune::VTK::Precision coordPrecision  = Dune::VTK::Precision::float64 ;
  // N.B. does not seem to be able to be allocated with new operator.
  Opm::DamarisVisOutput::DamarisVisDataWriter damarisGeomWriter(gridView, Dune::VTK::conforming, coordPrecision) ;
  // damarisGeomWriter = new DamarisVisDataWriter(gridView, Dune::VTK::conforming);  // this does not compile
  damarisGeomWriter.setupGeomData() ;
  damarisGeomWriter.printGridDetails() ;
 
   
*/

namespace Opm::DamarisVisOutput
{


  /**
   * @brief Writer for the ouput of grid functions in the vtk format.
   * @ingroup VTK
   *
   * Writes arbitrary grid functions (living on cells or vertices of a grid)
   * to a file suitable for easy visualization with
   * <a href="http://public.kitware.com/VTK/">The Visualization Toolkit (VTK)</a>.
   */
  template< class GridView >
  class DamarisVisDataWriter {



    // extract types
    typedef typename GridView::Grid Grid;
    typedef typename GridView::ctype DT;
    enum { n = GridView::dimension };
    enum { w = GridView::dimensionworld };

    typedef typename GridView::template Codim< 0 >::Entity Cell;
    typedef typename GridView::template Codim< n >::Entity Vertex;
    typedef Cell Entity;

    typedef typename GridView::IndexSet IndexSet;

    static const Dune::PartitionIteratorType VTK_Partition = Dune::InteriorBorder_Partition;
    //static const PartitionIteratorType VTK_Partition = All_Partition;

    typedef typename GridView::template Codim< 0 >
    ::template Partition< VTK_Partition >::Iterator
    GridCellIterator;
    typedef typename GridView::template Codim< n >
    ::template Partition< VTK_Partition >::Iterator
    GridVertexIterator;

    typedef typename GridCellIterator::Reference EntityReference;

    typedef typename GridView::template Codim< 0 >
    ::Entity::Geometry::LocalCoordinate Coordinate;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper< GridView > VertexMapper;

    // return true if entity should be skipped in Vertex and Corner iterator
    static bool skipEntity( const Dune::PartitionType entityType )
    {
      switch( VTK_Partition )
      {
        // for All_Partition no entity has to be skipped
        case Dune::All_Partition:             return false;
        case Dune::InteriorBorder_Partition:  return ( entityType != Dune::InteriorEntity );
        default: DUNE_THROW(Dune::NotImplemented,"Add check for this partition type");
      }
      return false ;
    }


  

  public:

    
    //! Iterator over the grids elements
    /**
     * This class iterates over the gridview's elements.  It is the same as
     * the gridview's Codim<0>::Iterator for the InteriorBorder_Partition,
     * except that it add a position() method.
     */
    class CellIterator : public GridCellIterator
    {
    public:
      //! construct a CellIterator from the gridview's Iterator.
      CellIterator(const GridCellIterator & x) : GridCellIterator(x) {}
      //! get the position of the center of the element, in element-local
      //! coordinates
      const Dune::FieldVector<DT,n> position() const
      {
        return Dune::Geo::ReferenceElements<DT,n>::general((*this)->type()).position(0,0);
      }
    };

    CellIterator cellBegin() const
    {
      return gridView_.template begin< 0, VTK_Partition >();
    }

    CellIterator cellEnd() const
    {
      return gridView_.template end< 0, VTK_Partition >();
    }

    //! Iterate over the grid's vertices
    /**
     * This class iterates over the elements, and within the elements over the
     * corners.  If the data mode dm is nonconforming, each vertex is visited
     * once for each element where it is a corner (similar to CornerIterator).
     * If dm is conforming each vertex is visited only once globally, for the
     * first element where it is a corner.  Contrary to CornerIterator, visit
     * the corners of a given element in Dune-ordering.
     *
     * Dereferencing the iterator yields the current entity, and the index of
     * the current corner within that entity is returned by the iterators
     * localindex() method.  Another useful method on the iterator itself is
     * position() which returns the element-local position of the current
     * corner.
     */
    class VertexIterator :
      public Dune::ForwardIteratorFacade<VertexIterator, const Entity, EntityReference, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      Dune::VTK::DataMode datamode;
      // Index of the currently visited corner within the current element.
      // NOTE: this is in Dune-numbering, in contrast to CornerIterator.
      int cornerIndexDune;
      const VertexMapper & vertexmapper;
      std::vector<bool> visited;
      // in conforming mode, for each vertex id (as obtained by vertexmapper)
      // hold its number in the iteration order (VertexIterator)
      int offset;

      // hide operator ->
      void operator->();
    protected:
      void basicIncrement ()
      {
        if( git == gend )
          return;
        ++cornerIndexDune;
        const int numCorners = git->subEntities(n);
        if( cornerIndexDune == numCorners )
        {
          offset += numCorners;
          cornerIndexDune = 0;

          ++git;
          while( (git != gend) && skipEntity( git->partitionType() ) )
            ++git;
        }
      }
    public:
      VertexIterator(const GridCellIterator & x,
                     const GridCellIterator & end,
                     const Dune::VTK::DataMode & dm,
                     const VertexMapper & vm) :
        git(x), gend(end), datamode(dm), cornerIndexDune(0),
        vertexmapper(vm), visited(vm.size(), false),
        offset(0)
      {
        if (datamode == Dune::VTK::conforming && git != gend)
          visited[vertexmapper.subIndex(*git,cornerIndexDune,n)] = true;
      }
      void increment ()
      {
        switch (datamode)
        {
        case Dune::VTK::conforming :
          while(visited[vertexmapper.subIndex(*git,cornerIndexDune,n)])
          {
            basicIncrement();
            if (git == gend) return;
          }
          visited[vertexmapper.subIndex(*git,cornerIndexDune,n)] = true;
          break;
        case Dune::VTK::nonconforming :
          basicIncrement();
          break;
        }
      }
      bool equals (const VertexIterator & cit) const
      {
        return git == cit.git
               && cornerIndexDune == cit.cornerIndexDune
               && datamode == cit.datamode;
      }
      EntityReference dereference() const
      {
        return *git;
      }
      //! index of vertex within the entity, in Dune-numbering
      int localindex () const
      {
        return cornerIndexDune;
      }
      //! position of vertex inside the entity
      Dune::FieldVector<DT,n> position () const
      {
        return Dune::referenceElement<DT,n>(git->type())
          .position(cornerIndexDune,n);
      }
    };

    VertexIterator vertexBegin () const
    {
      return VertexIterator( gridView_.template begin< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper );
    }

    VertexIterator vertexEnd () const
    {
      return VertexIterator( gridView_.template end< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper );
    }

    //! Iterate over the elements' corners
    /**
     * This class iterates over the elements, and within the elements over the
     * corners.  Each vertex in the grid can be a corner in multiple elements,
     * and is visited once for each element it is associated with.  This class
     * differs from VertexIterator in that it visits the corners of a given
     * element in VTK-ordering, and that it always visits a given vertex once
     * for each element where that vertex is a corner in, independent of the
     * data mode dm.
     *
     * Dereferencing the iterator yields the current entity.  Another useful
     * method on the iterator itself is id(), which returns the number of the
     * current corners associated vertex, in the numbering given by the
     * iteration order of VertexIterator.
     */
    class CornerIterator :
      public Dune::ForwardIteratorFacade<CornerIterator, const Entity, EntityReference, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      Dune::VTK::DataMode datamode;
      // Index of the currently visited corner within the current element.
      // NOTE: this is in VTK-numbering, in contrast to VertexIterator.
      int cornerIndexVTK;
      const VertexMapper & vertexmapper;
      // in conforming mode, for each vertex id (as obtained by vertexmapper)
      // hold its number in the iteration order of VertexIterator (*not*
      // CornerIterator)
      const std::vector<int> & number;
      // holds the number of corners of all the elements we have seen so far,
      // excluding the current element
      int offset;

      // hide operator ->
      void operator->();
    public:
      CornerIterator(const GridCellIterator & x,
                     const GridCellIterator & end,
                     const Dune::VTK::DataMode & dm,
                     const VertexMapper & vm,
                     const std::vector<int> & num) :
        git(x), gend(end), datamode(dm), cornerIndexVTK(0),
        vertexmapper(vm),
        number(num), offset(0) {}
      void increment ()
      {
        if( git == gend )
          return;
        ++cornerIndexVTK;
        const int numCorners = git->subEntities(n);
        if( cornerIndexVTK == numCorners )
        {
          offset += numCorners;
          cornerIndexVTK = 0;

          ++git;
          while( (git != gend) && skipEntity( git->partitionType() ) )
            ++git;
        }
      }
      bool equals (const CornerIterator & cit) const
      {
        return git == cit.git
               && cornerIndexVTK == cit.cornerIndexVTK
               && datamode == cit.datamode;
      }
      EntityReference dereference() const
      {
        return *git;
      }
      //! Process-local consecutive zero-starting vertex id
      /**
       * This method returns the number of this corners associated vertex, in
       * the numbering given by the iteration order of VertexIterator.
       */
      int id () const
      {
        switch (datamode)
        {
        case Dune::VTK::conforming :
          return
            number[vertexmapper.subIndex(*git,Dune::VTK::renumber(*git,cornerIndexVTK),
                                    n)];
        case Dune::VTK::nonconforming :
          return offset + Dune::VTK::renumber(*git,cornerIndexVTK);
        default :
          DUNE_THROW(Dune::IOError,"DamarisVisDataWriter: unsupported DataMode" << datamode);
        }
      }
    };

    CornerIterator cornerBegin () const
    {
      return CornerIterator( gridView_.template begin< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }

    CornerIterator cornerEnd () const
    {
      return CornerIterator( gridView_.template end< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }

  public:
    /**
     * @brief Construct a DamarisVisDataWriter working on a specific GridView.
     *
     *
     * @param gridView The gridView the grid functions live on. (E. g. a LevelGridView.)
     * @param dm The data mode.
     * @param coordPrecision the precision with which to write out the coordinates
     */
    explicit DamarisVisDataWriter ( const GridView &gridView,
                         Dune::VTK::DataMode dm = Dune::VTK::conforming)
      : gridView_( gridView ),
        datamode( dm ),
       /* coordPrec (coordPrecision),*/
        polyhedralCellsPresent_( checkForPolyhedralCells() )
    { 
        // Storing the expected name of the variables in the Damaris XML file
        coordset_coords_values_x = "coordset/coords/values/x";
        coordset_coords_values_y = "coordset/coords/values/y";
        coordset_coords_values_z = "coordset/coords/values/z";
        
        topologies_topo_elements_connectivity = "topologies/topo/elements/connectivity";
        topologies_topo_elements_offsets      = "topologies/topo/elements/offsets";
        topologies_topo_elements_types        = "topologies/topo/elements/types";
        
        // Only needed if polyhedralCellsPresent_ == true
        topologies_topo_subelements_faces = "topologies/topo/subelements/faces";
        topologies_topo_subelements_offsets  = "topologies/topo/subelements/offsets";
        
        this->setup_called_bool_ = false ;
        this->setupGeomData() ;
    
    }


    //! clear list of registered functions
    void clear ()
    {
      //celldata.clear();
      //vertexdata.clear();
    }

    //! get the precision with which coordinates are written out
    /*Dune::VTK::Precision coordPrecision() const
    { return coordPrec; }
*/
    //! destructor
    ~DamarisVisDataWriter ()
    {
      this->clear();
    }

    

    //! Get the detais of the grid data 
    void setupGeomData ( void )
    {
      Dune::VTK::FileType fileType =
        (n == 1) ? Dune::VTK::polyData : Dune::VTK::unstructuredGrid;

      // VTK::VTUWriter writer(s, outputtype, fileType);

      // Grid characteristics
      vertexmapper = new VertexMapper( gridView_, Dune::mcmgVertexLayout() );
      if (datamode == Dune::VTK::conforming)
      {
        number.resize(vertexmapper->size());
        for (std::vector<int>::size_type i=0; i<number.size(); i++) 
            number[i] = -1;
      }
      countEntities(nvertices, ncells, ncorners);

      this->setup_called_bool_ = true ;

      delete vertexmapper; number.clear();
    }

    void writeGeometryData() {
      // PointData fields
      //writeVertexData(writer);

      // CellData fields
      //writeCellData(writer);

      // x,y,z vertices coordinates
      // writeGridPoints(writer);

      // Cells -connectivity, offested, types, polygon extras if required
      // writeGridCells(writer);
    }


    std::string getTypeString() const
    {
      if (n==1)
        return "PolyData";
      else
        return "UnstructuredGrid";
    }

    //! count the vertices, cells and corners
    void countEntities(int &nvertices_, int &ncells_, int &ncorners_)
    {
      nvertices_ = 0;
      ncells_ = 0;
      ncorners_ = 0;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        ncells_++;
        // because of the use of vertexmapper->map(), this iteration must be
        // in the order of Dune's numbering.
        const int subEntities = it->subEntities(n);
        for (int i=0; i<subEntities; ++i)
        {
          ncorners_++;
          if (datamode == Dune::VTK::conforming)
          {
            int alpha = vertexmapper->subIndex(*it,i,n);
            if (number[alpha]<0)
              number[alpha] = nvertices_++;
          }
          else
          {
            nvertices_++;
          }
        }
      }
    }

/*
    template<typename T>
    std::tuple<std::string,std::string> getDataNames(const T& data) const
    {
      std::string scalars = "";
      for (auto it = data.begin(),
             end = data.end();
           it != end;
           ++it)
        if (it->fieldInfo().type() == VTK::FieldInfo::Type::scalar)
          {
            scalars = it->name();
            break;
          }

      std::string vectors = "";
      for (auto it = data.begin(),
             end = data.end();
           it != end;
           ++it)
        if (it->fieldInfo().type() == VTK::FieldInfo::Type::vector)
          {
            vectors = it->name();
            break;
          }
      return std::make_tuple(scalars,vectors);
    }

    template<typename Data, typename Iterator>
    void writeData(VTK::VTUWriter& writer, const Data& data, const Iterator begin, const Iterator end, int nentries)
    {
      for (auto it = data.begin(),
             iend = data.end();
           it != iend;
           ++it)
      {
        const auto& f = *it;
        VTK::FieldInfo fieldInfo = f.fieldInfo();
        std::size_t writecomps = fieldInfo.size();
        switch (fieldInfo.type())
          {
          case VTK::FieldInfo::Type::scalar:
            break;
          case VTK::FieldInfo::Type::vector:
            // vtk file format: a vector data always should have 3 comps (with
            // 3rd comp = 0 in 2D case)
            if (writecomps > 3)
              DUNE_THROW(IOError,"Cannot write VTK vectors with more than 3 components (components was " << writecomps << ")");
            writecomps = 3;
            break;
          case VTK::FieldInfo::Type::tensor:
            DUNE_THROW(NotImplemented,"VTK output for tensors not implemented yet");
          }
        std::shared_ptr<VTK::DataArrayWriter> p
          (writer.makeArrayWriter(f.name(), writecomps, nentries, fieldInfo.precision()));
        if(!p->writeIsNoop())
          for (Iterator eit = begin; eit!=end; ++eit)
          {
            const Entity & e = *eit;
            f.bind(e);
            f.write(eit.position(),*p);
            f.unbind();
            // vtk file format: a vector data always should have 3 comps
            // (with 3rd comp = 0 in 2D case)
            for (std::size_t j=fieldInfo.size(); j < writecomps; ++j)
              p->write(0.0);
          }
      }
    }

    //! write cell data
    virtual void writeCellData(VTK::VTUWriter& writer)
    {
      if(celldata.size() == 0)
        return;

      std::string scalars, vectors;
      std::tie(scalars,vectors) = getDataNames(celldata);

      writer.beginCellData(scalars, vectors);
      writeData(writer,celldata,cellBegin(),cellEnd(),ncells);
      writer.endCellData();
    }

    //! write vertex data
    virtual void writeVertexData(VTK::VTUWriter& writer)
    {
      if(vertexdata.size() == 0)
        return;

      std::string scalars, vectors;
      std::tie(scalars,vectors) = getDataNames(vertexdata);

      writer.beginPointData(scalars, vectors);
      writeData(writer,vertexdata,vertexBegin(),vertexEnd(),nvertices);
      writer.endPointData();
    }*/

    //! write the positions of vertices
    template <typename T>
    void writeGridPoints( std::vector<T> & x_inout,  std::vector<T> & y_inout, std::vector<T> & z_inout )
    {
        int dimw=w;
        VertexIterator vEnd = vertexEnd();
        for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit)
        {          
          x_inout.push_back(  (T) (*vit).geometry().corner(vit.localindex())[0] );
          y_inout.push_back(  (T) (*vit).geometry().corner(vit.localindex())[1] );
          if (dimw == 3) 
            z_inout.push_back( (T) (*vit).geometry().corner(vit.localindex())[2] ); 
          else 
            z_inout.push_back((T)  0.0);  
        }
    }
    
    //! write the positions of vertices
    template <typename T>
    void writeGridPoints( T*  x_inout,  T*  y_inout, T* z_inout )
    {
        int dimw=w;
        VertexIterator vEnd = vertexEnd();
        int i = 0 ;
        T zero_val = 0.0 ;
        for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit)
        {          
          x_inout[i] = static_cast<T>((*vit).geometry().corner(vit.localindex())[0]) ;
          y_inout[i] = static_cast<T>((*vit).geometry().corner(vit.localindex())[1]) ;
          if (dimw == 3) 
            z_inout[i] = static_cast<T>((*vit).geometry().corner(vit.localindex())[2]) ; 
          else 
            z_inout[i] = zero_val;  
        
          i++ ;
        }
    }

    //! write the connectivity array
    template <typename T, typename U, typename V, typename W>
    void writeGridCells(std::vector<T> & connectivity_inout,  std::vector<U> & offsets_inout, std::vector<V> & types_inout, std::vector<W> & polyfaces_inout, std::vector<W> & polyfacesoffsets_inout  )
    {
      // connectivity
       for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
              connectivity_inout.push_back(static_cast<T>(it.id())) ;

      // offsets
      U offset = 0;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        offset += static_cast<U>(it->subEntities(n));
        offsets_inout.push_back(offset) ;
      }


      // types
      if (n>1)
      {
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          int vtktype = Dune::VTK::geometryType(it->type());
            types_inout.push_back(static_cast<V>(vtktype));
        }
        
        // if polyhedron cells found also cell faces need to be written
        if( polyhedralCellsPresent_ )
        {
          writeCellFaces( polyfaces_inout, polyfacesoffsets_inout );
        }
      }
    }
    
        /**
    * write the connectivity array
    */
    template <typename T>
    void writeConnectivity(std::vector<T> & connectivity_inout )
    {
    
      // connectivity
       int i = 0 ;
       for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
       {
           T connect_data = static_cast<T>(it.id()) ;
           connectivity_inout.push_back( connect_data );
       }

    }
    
    /**
    * write the offsets values
    */
    template <typename T>
    void writeOffsetsCells(std::vector<T> & offsets_inout  )
    {
      // offsets
      T offset = 0;
      int i = 0 ;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        offset += static_cast<T>(it->subEntities(n));
        offsets_inout.push_back( offset ) ;
      }
    }
    
    /**
    * write the Cell types array
    */
    template <typename T>
    void writeCellTypes( std::vector<T> & types_inout )
    {
      int i = 0 ;
      // types
      if (n>1)
      {
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          T vtktype = static_cast<T>(Dune::VTK::geometryType(it->type()));
          types_inout.push_back( vtktype ) ;
        }
      }
    }
    

    /**
    * write the connectivity array
    */
    template <typename T>
    void writeConnectivity(T * connectivity_inout )
    {
    
      // connectivity
       int i = 0 ;
       for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
       {
           T connect_data = static_cast<T>(it.id()) ;
           connectivity_inout[i++] = connect_data ;
       }

    }
    
    /**
    * write the offsets values
    */
    template <typename T>
    void writeOffsetsCells( T* offsets_inout  )
    {
      // offsets
      T offset = 0;
      int i = 0 ;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        offset += static_cast<T>(it->subEntities(n));
        offsets_inout[i++] = offset ;
      }
    }
    
    /**
    * write the Cell types array
    */
    template <typename T>
    void writeCellTypes( T* types_inout )
    {
      int i = 0 ;
      // types
      if (n>1)
      {
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          T vtktype = static_cast<T>(Dune::VTK::geometryType(it->type()));
          types_inout[i++] = vtktype ;

        }
      }
    }
    
    
        /**
    * write the connectivity array
    * N.B. W* polyfaces_inout, X* polyfacesoffsets_inout are only needed if hasPolyhedralCells() returns true
    */
    template <typename T>
    void writeGridCells( T* polyfaces_inout, T* polyfacesoffsets_inout  )
    {
    
      int i = 0 ;
      // types
      if (n>1)
      {
        // if polyhedron cells found also cell faces need to be written
        if( polyhedralCellsPresent_ )
        {
          writeCellFaces( polyfaces_inout, polyfacesoffsets_inout );
        }
      }

    }
    

    bool checkForPolyhedralCells() const
    {
      // check if polyhedron cells are present
      for( const auto& geomType : gridView_.indexSet().types( 0 ) )
      {
        if( Dune::VTK::geometryType( geomType ) == Dune::VTK::polyhedron )
        {
          return true;
        }
      }
      return false;
    }

    //! write the connectivity array
    template <typename T>
    void writeCellFaces( std::vector<T> & polyfaces_inout, std::vector<T> & polyfacesoffsets_inout )
    {
      if( ! faceVertices_ )
      {
        faceVertices_.reset( new std::pair< std::vector<T>, std::vector<T> > () );
        // fill face vertex structure
        fillFaceVertices( cornerBegin(), cornerEnd(), gridView_.indexSet(),
                          faceVertices_->first, faceVertices_->second );
      }

      std::vector< T >& faces = faceVertices_->first;
      std::vector< T >& faceOffsets = faceVertices_->second;
      assert( int(faceOffsets.size()) == ncells );

      {
          T face_typew ;
          for( const auto& face : faces )
          {
              face_typew = static_cast<T>(face);
              polyfaces_inout.push_back(face) ;
          }
      }

      {
          T faceoffset_typeX ;
          for( const auto& offset : faceOffsets )
          {
              faceoffset_typeX = static_cast<T>(offset);
              polyfacesoffsets_inout.push_back(faceoffset_typeX) ;
          }
          // clear face vertex structure
          faceVertices_.reset();
      }
    }

    template <class CornerIterator, class IndexSet, class T>
    inline void fillFaceVertices( CornerIterator it,
                           const CornerIterator end,
                           const IndexSet& indexSet,
                           std::vector<T>& faces,
                           std::vector<T>& faceOffsets )
    {
      if( n == 3 && it != end )
      {
        // clear output arrays
        faces.clear();
        faces.reserve( 15 * ncells );
        faceOffsets.clear();
        faceOffsets.reserve( ncells );

        int offset = 0;

        Cell element = *it;
        int elIndex = indexSet.index( element );
        std::vector< T > vertices;
        vertices.reserve( 30 );
        for( ; it != end; ++it )
        {
          const Cell& cell = *it ;
          const int cellIndex = indexSet.index( cell ) ;
          if( elIndex != cellIndex )
          {
            fillFacesForElement( element, indexSet, vertices, offset, faces, faceOffsets );

            vertices.clear();
            element = cell ;
            elIndex = cellIndex ;
          }
          vertices.push_back( it.id() );
        }

        // fill faces for last element
        fillFacesForElement( element, indexSet, vertices, offset, faces, faceOffsets );
      }
    }

    template <class Entity, class IndexSet, class T>
    static void fillFacesForElement( const Entity& element,
                                     const IndexSet& indexSet,
                                     const std::vector<T>& vertices,
                                     T& offset,
                                     std::vector<T>& faces,
                                     std::vector<T>& faceOffsets )
    {
      const int dim = n;

      std::map< T, T > vxMap;

      // get number of local faces
      const int nVertices = element.subEntities( dim );
      for( int vx = 0; vx < nVertices; ++ vx )
      {
        const int vxIdx = indexSet.subIndex( element, vx, dim );
        vxMap[ vxIdx ] = vertices[ vx ];
      }

      // get number of local faces
      const int nFaces = element.subEntities( 1 );
      // store number of faces for current element
      faces.push_back( nFaces );
      ++offset;
      // extract each face as a set of vertex indices
      for( int fce = 0; fce < nFaces; ++ fce )
      {
        // obtain face
        const auto face = element.template subEntity< 1 > ( fce );

        // get all vertex indices from current face
        const int nVxFace = face.subEntities( dim );
        faces.push_back( nVxFace );
        ++offset ;
        for( int i=0; i<nVxFace; ++i )
        {
          const T vxIndex = indexSet.subIndex( face, i, dim );
          assert( vxMap.find( vxIndex ) != vxMap.end() );
          faces.push_back( vxMap[ vxIndex ] );
          ++offset ;
        }
      }

      // store face offset for each element
      faceOffsets.push_back( offset );
    }
    
  void   printGridDetails()
  {
      printNCells() ;
      printNVertices() ;
      printNCorners() ;
      std::cout << "Mesh Type = " << getTypeString()  << std::endl ;
  }
    
  void printNCells()
  {
      std::cout << "ncells = " << ncells << std::endl ;
  }
  
  void printNVertices()
  {
      std::cout << "nvertices = " << nvertices << std::endl ;
  }
  
  void printNCorners()
  {
      std::cout << "ncorners = " << ncorners << std::endl ;
  }
  
  int getNCells()
  {
      return(ncells) ;
  }
  
  int getNVertices()
  {
      return(nvertices) ;
  }
  
  int getNCorners()
  {
      return(ncorners) ;
  }
  
  bool hasPolyhedralCells( void )
  {
      return (polyhedralCellsPresent_) ;
  }
  
  /** 
  * Used if we need to change the variable name as used in the Damaris XML file
  * The defaults are set in the DamarisVisDataWriter constructor.
  *
  *  Vertex data
  */
  void set_damaris_var_name_vertex_x( std::string& var_name )
  {
      coordset_coords_values_x =  var_name ;
  }
  void set_damaris_var_name_vertex_y( std::string& var_name )
  {
      coordset_coords_values_y =  var_name ;
  }
  void set_damaris_var_name_vertex_z( std::string& var_name )
  {
      coordset_coords_values_z =  var_name ;
  }
  
  /** 
  * Used if we need to change the variable name as used in the Damaris XML file
  * The defaults are set in the DamarisVisDataWriter constructor.
  *
  *  Connectivity data
  */
  void set_damaris_var_name_connectivity( std::string& var_name )
  {
      topologies_topo_elements_connectivity =  var_name ;
  }
  void set_damaris_var_name_offsets( std::string& var_name )
  {
      topologies_topo_elements_offsets =  var_name ;
  }
  void set_damaris_var_name_types( std::string& var_name )
  {
      topologies_topo_elements_types =  var_name ;
  }
  
  /** 
  * Used if we need to change the variable name as used in the Damaris XML file
  * The defaults are set in the DamarisVisDataWriter constructor.
  *
  *  Polyhedral Data
  */ 
  void set_damaris_var_name_poly_faces( std::string& var_name )
  {
      topologies_topo_subelements_faces =  var_name ;
  }
  void set_damaris_var_name_poly_offsets( std::string& var_name )
  {
      topologies_topo_subelements_offsets =  var_name ;
  }
  
   /** 
  * Used if we need to change the variable name as used in the Damaris XML file
  * The defaults are set in the DamarisVisDataWriter constructor.
  *
  *  Vertex data
  */
  void printDamarisVarNames( int rank)
  {
      std::cout << "Rank:" << rank << " coordset_coords_values_x = " << coordset_coords_values_x << std::endl ; ;

      std::cout << "Rank:" << rank << " coordset_coords_values_y  = " << coordset_coords_values_y << std::endl ;

      std::cout << "Rank:" << rank << " coordset_coords_values_z = " << coordset_coords_values_z << std::endl ;

      std::cout << "Rank:" << rank << " topologies_topo_elements_connectivity = " << topologies_topo_elements_connectivity << std::endl ;

      std::cout << "Rank:" << rank << " topologies_topo_elements_offsets = " << topologies_topo_elements_offsets << std::endl ;

      std::cout << "Rank:" << rank << " topologies_topo_elements_types = " << topologies_topo_elements_types << std::endl ;

      std::cout << "Rank:" << rank << " topologies_topo_subelements_faces = " << topologies_topo_subelements_faces << std::endl ;

      std::cout << "Rank:" << rank << " topologies_topo_subelements_offsets = " << topologies_topo_subelements_offsets << std::endl ;
  }
  
  template <typename T>
  void SetPointersToDamarisShmem( std::string damaris_variable_name, T** ret_ptr, int rank )
  {
     int nverticies = this->getNVertices();
     int damaris_err ;

    T * temp_ptr ;
    // Allocate memory in the shared memory section... 
    damaris_err = damaris_alloc(damaris_variable_name.c_str(), (void **) &temp_ptr) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::SetDataVertexPointers: damaris_alloc(\"" << damaris_variable_name <<"\", (void **) &ret_ptr)" << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    *ret_ptr = temp_ptr ;
  }
  
  
  
  void SetDamarisParameter(std::string param_name, int paramSizeVal, int rank) 
  {
    int damaris_err ;
    damaris_err = damaris_parameter_set(param_name.c_str(), &paramSizeVal, sizeof(int));
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::SetDataVertexPointers : damaris_parameter_set(\"" << param_name << "\", paramSizeVal, sizeof(int));  Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
  }
  
  template <typename T>
  void SetDataVertexPointers( T ** x_dptr, T ** y_dptr, T ** z_dptr, int rank )
  {
    
    int nverticies = this->getNVertices();
    // this sets the size of the shared memory block obtained
    SetDamarisParameter("n_coords_local", nverticies, rank);
    
    SetPointersToDamarisShmem(coordset_coords_values_x, x_dptr, rank) ;
    
    SetPointersToDamarisShmem(coordset_coords_values_y, y_dptr, rank) ;
    
    SetPointersToDamarisShmem(coordset_coords_values_z, z_dptr, rank) ;
        
  }
  
  void VertexArraysCommitAndClear( int rank ) 
  {
    int damaris_err ;
    // Signal to Damaris we are done writing data for this iteration
    damaris_err = damaris_commit (coordset_coords_values_x.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_commit(\"" << coordset_coords_values_x <<"\")"  << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    damaris_err = damaris_commit (coordset_coords_values_y.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_commit(\"" << coordset_coords_values_y <<"\")"  << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    damaris_err = damaris_commit (coordset_coords_values_z.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_commit(\"" << coordset_coords_values_z <<"\")"  << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    
    // Signal to Damaris we are done accessing data for this iteration
    damaris_err = damaris_clear(coordset_coords_values_x.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_clear(\"" << coordset_coords_values_x << "\")"  << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    damaris_err = damaris_clear(coordset_coords_values_y.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_clear(\"" << coordset_coords_values_y << "\")" <<  ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    damaris_err = damaris_clear(coordset_coords_values_z.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_clear(\"" << coordset_coords_values_z << "\")" << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }  
  }
  

  void ConnectivityArraysCommitAndClear( int rank ){
    int damaris_err ;
    // Signal to Damaris we are done writing data for this iteration
    damaris_err = damaris_commit (topologies_topo_elements_connectivity.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_commit(\"" << topologies_topo_elements_connectivity <<"\")"  << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    damaris_err = damaris_commit (topologies_topo_elements_offsets.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_commit(\"" << topologies_topo_elements_offsets <<"\")"  << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    damaris_err = damaris_commit (topologies_topo_elements_types.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_commit(\"" << topologies_topo_elements_types <<"\")"  << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    
    // Signal to Damaris we are done accessing data for this iteration
    damaris_err = damaris_clear(topologies_topo_elements_connectivity.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_clear(\"" << topologies_topo_elements_connectivity << "\")"  << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    damaris_err = damaris_clear(topologies_topo_elements_offsets.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_clear(\"" << topologies_topo_elements_offsets << "\")" <<  ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }
    damaris_err = damaris_clear(topologies_topo_elements_types.c_str()) ;
    if (damaris_err != DAMARIS_OK) {
        std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::VertexArraysCommitAndClear : damaris_clear(\"" << topologies_topo_elements_types << "\")" << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
    }  
  }
  
  void PassVertexDataToDamaris( int rank )
  {
    // We assume setupGeomData() has been called
    if (this->setup_called_bool_ == true) {
            
        int damaris_err ;
        // Get the expected data type of the vertex data 
        // This is defined in the Damaris XML file layout for the vertex variables.
        // We are assuming the Damaris vertex variables (y and z coordinates) will
        // have the same type as the x vertex data.
        DAMARIS_TYPE_STR vartype ;
        damaris_err = damaris_get_type(coordset_coords_values_x.c_str(), &vartype)  ;
        if (damaris_err != DAMARIS_OK) {
            std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::PassVertexDataToDamaris : damaris_alloc(\"coordsets/coords/values/x\", (void **) &x_dptr)" << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
        }

        if (vartype == DAMARIS_TYPE_DOUBLE)
        {
            std::cout << "INFO: The damaris variable has type \"double\"" << std::endl ;
            double * x_dptr ;
            double * y_dptr ;
            double * z_dptr ;
            
            this->SetDataVertexPointers(&x_dptr, &y_dptr, &z_dptr, rank) ;

            // Write the data to the Damaris shared memory
            this->writeGridPoints(x_dptr, y_dptr, z_dptr) ;
            
        }
        else if (vartype == DAMARIS_TYPE_FLOAT)
        {
            std::cout << "INFO: The damaris vertex variable has type \"float\"" << std::endl ;
            float * x_fptr ;
            float * y_fptr ;
            float * z_fptr ;
            
            // This allocates space in the Damaris shared memory region
            this->SetDataVertexPointers(&x_fptr, &y_fptr, &z_fptr, rank) ;

            // Write the data to the Damaris shared memory
            this->writeGridPoints(x_fptr, y_fptr, z_fptr) ;
        }
        // Tell Damaris that we are finished writing (or reading) the vertex data
        this->VertexArraysCommitAndClear( rank ) ;

    } else {
        std::cerr << "ERROR = " << rank << " : DamarisVisDataWriter::PassVertexDataToDamaris(): setupGeomData() has not been called" << std::endl ;
        // exit(-1) ;
    }
  }
  
  
  void PassConnectivityDataToDamaris( int rank )
  {
    // We assume setupGeomData() has been called
    if (this->setup_called_bool_ == true) {
            
        int damaris_err ;
        // Get the expected data type of the vertex data 
        // This is defined in the Damaris XML file layout for the vertex variables.
        // We are assuming the Damaris vertex variables (y and z coordinates) will
        // have the same type as the x vertex data.
        DAMARIS_TYPE_STR connectvartype ;
        damaris_err = damaris_get_type(topologies_topo_elements_connectivity.c_str(), &connectvartype)  ;
        if (damaris_err != DAMARIS_OK) {
            std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::PassConnectivityDataToDamaris : damaris_get_type(\"" << topologies_topo_elements_connectivity << "\", connectvartype)" << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
        }

        if (connectvartype == DAMARIS_TYPE_INT)
        {
            std::cout << "INFO: The damaris connectvartype has type \"int\"" << std::endl ;
            int * connect_intptr ;
            int * offset_intptr ;
            std::vector<int> connect_vect ;
            std::vector<int> offset_vect ;
            // this->SetDataConnectionPointers(&connect_intptr) ;
            
            this->writeConnectivity(connect_vect) ;
            int paramSizeVal = connect_vect.size() ;
            this->SetDamarisParameter("n_connectivity_ph",  paramSizeVal,  rank) ;
            // This allocates space in the Damaris shared memory region
            this->SetPointersToDamarisShmem(topologies_topo_elements_connectivity, &connect_intptr, rank) ;    
            std::cout << "Rank : " << rank << " INFO: The damaris connect_vect array has size = " << connect_vect.size() << std::endl ;            
            for (int i = 0 ; i < connect_vect.size() ; i++)
            {
                connect_intptr[i] = connect_vect[i] ;
            }
            
            this->writeOffsetsCells(offset_vect) ;
            paramSizeVal = offset_vect.size() ;
            this->SetDamarisParameter("n_offsets_types_ph",  paramSizeVal,  rank) ;
            // This allocates space in the Damaris shared memory region
            this->SetPointersToDamarisShmem(topologies_topo_elements_offsets, &offset_intptr, rank) ;
            std::cout << "Rank : " << rank << " INFO: The damaris offset_vect array has size = " << offset_vect.size() << std::endl ;      
            for (int i = 0 ; i < offset_vect.size() ; i++)
            {
                offset_intptr[i] = offset_vect[i] ;
            }
            
        }
        else if (connectvartype == DAMARIS_TYPE_LONG)
        {
            std::cout << "INFO: The damaris connectvartype variable has type \"long\"" << std::endl ;
            long * connect_longptr ;
            long * offset_longptr ;
            std::vector<long> connect_vect ;
            std::vector<long> offset_vect ;
            // this->SetDataConnectionPointers(&connect_intptr) ;
            
            this->writeConnectivity(connect_vect) ;
            int paramSizeVal = connect_vect.size() ;
            this->SetDamarisParameter("n_connectivity_ph",  paramSizeVal,  rank) ;
            // This allocates space in the Damaris shared memory region
            this->SetPointersToDamarisShmem(topologies_topo_elements_connectivity, &connect_longptr, rank) ;   
            std::cout << "Rank : " << rank << " INFO: The damaris connect_vect array has size = " << connect_vect.size() << std::endl ;                  
            for (int i = 0 ; i < connect_vect.size() ; i++)
            {
                connect_longptr[i] = connect_vect[i] ;
            }
            
            this->writeOffsetsCells(offset_vect) ;
            paramSizeVal = offset_vect.size() ;
            this->SetDamarisParameter("n_offsets_types_ph",  paramSizeVal,  rank) ;
            // This allocates space in the Damaris shared memory region
            this->SetPointersToDamarisShmem(topologies_topo_elements_offsets, &offset_longptr, rank) ;
            std::cout << "Rank : " << rank << " INFO: The damaris offset_vect array has size = " << offset_vect.size() << std::endl ;      
            for (int i = 0 ; i < offset_vect.size() ; i++)
            {
                offset_longptr[i] = offset_vect[i] ;
            }
        }
        
        
        // Set up the cell types arrays and pass to Damaris
        DAMARIS_TYPE_STR celltypevartype ;
        damaris_err = damaris_get_type(topologies_topo_elements_types.c_str(), &celltypevartype)  ;
        if (damaris_err != DAMARIS_OK) {
            std::cerr << "ERROR rank =" << rank << " : DamarisVisDataWriter::PassConnectivityDataToDamaris : damaris_get_type(\"" << topologies_topo_elements_types << "\", celltypevartype)" << ", Damaris error = " <<  damaris_error_string(damaris_err) << std::endl ;
        }

        if (celltypevartype == DAMARIS_TYPE_CHAR)
        {
            std::cout << "INFO: The damaris celltypevartype has type \"char\"" << std::endl ;
            char * celltype_charptr ;
            std::vector<char> celltype_vect ;
            this->writeCellTypes(celltype_vect) ;
            int paramSizeVal = celltype_vect.size() ;
            this->SetDamarisParameter("n_offsets_types_ph",  paramSizeVal,  rank) ;  // this is re-set as it is set above
            // This allocates space in the Damaris shared memory region
            this->SetPointersToDamarisShmem(topologies_topo_elements_types, &celltype_charptr, rank) ;    
            std::cout << "Rank : " << rank << " INFO: The damaris celltype_vect array has size = " << celltype_vect.size() << std::endl ; 
            for (int i = 0 ; i < celltype_vect.size() ; i++)
            {
                celltype_charptr[i] = celltype_vect[i] ;
            }
            
            
        }
        else if (celltypevartype == DAMARIS_TYPE_UCHAR)
        {
            std::cout << "INFO: The damaris celltypevartype has type \"unsigned char\"" << std::endl ;
            unsigned char * celltype_ucharptr ;
            std::vector<char> celltype_vect ;
            this->writeCellTypes(celltype_vect) ;
            int paramSizeVal = celltype_vect.size() ;
            this->SetDamarisParameter("n_offsets_types_ph",  paramSizeVal,  rank) ; // this is re-set as it is set above
            // This allocates space in the Damaris shared memory region
            this->SetPointersToDamarisShmem(topologies_topo_elements_types, &celltype_ucharptr, rank) ;        
            std::cout << "Rank : " << rank << " INFO: The damaris celltype_vect array has size = " << celltype_vect.size() << std::endl ;
            for (int i = 0 ; i < celltype_vect.size() ; i++)
            {
                celltype_ucharptr[i] = celltype_vect[i] ;
            }
        } else {
            std::cout << "ERROR: DamarisVisDataWriter::PassConnectivityDataToDamaris The damaris celltypevartype did was not recognized!" << std::endl ;
        }
        
       
        
        // Tell Damaris that we are finished writing (or reading) the vertex data
        this->ConnectivityArraysCommitAndClear( rank ) ;

    } else {
        std::cerr << "ERROR = " << rank << " : DamarisVisDataWriter::PassConnectivityDataToDamaris(): setupGeomData() has not been called" << std::endl ;
        // exit(-1) ;
    }
  }
  
  protected:
    // the list of registered functions
    // std::list<VTKLocalFunction> celldata;
    // std::list<VTKLocalFunction> vertexdata;

    // the grid
    GridView gridView_;

    // temporary grid information
    int ncells;
    int nvertices;
    int ncorners;
    
    // Storing the expected name of the variables in the Damaris XML file
    // These names will be used for introspection on the data type that Damaris 
    // is expecting which is defined in the Layouts of the variables in the XML file
    std::string coordset_coords_values_x ;
    std::string coordset_coords_values_y ;
    std::string coordset_coords_values_z ;
    
    std::string topologies_topo_elements_connectivity ;
    std::string topologies_topo_elements_offsets ;
    std::string topologies_topo_elements_types ;
    
    // Only needed if polyhedralCellsPresent_ == true
    std::string topologies_topo_subelements_faces ;
    std::string topologies_topo_subelements_offsets ;
  private:
    VertexMapper* vertexmapper;
    // in conforming mode, for each vertex id (as obtained by vertexmapper)
    // hold its number in the iteration order (VertexIterator)
    std::vector<int> number;
    Dune::VTK::DataMode datamode;
    // Dune::VTK::Precision coordPrec;

    // true if polyhedral cells are present in the grid
    const bool polyhedralCellsPresent_;
    
    // Used to check that class data has been initialized
    bool setup_called_bool_ ;

    // pointer holding face vertex connectivity if needed
    std::shared_ptr< std::pair< std::vector<int>, std::vector<int> > > faceVertices_;

  protected:
    Dune::VTK::OutputType outputtype;
  };

}

#endif
