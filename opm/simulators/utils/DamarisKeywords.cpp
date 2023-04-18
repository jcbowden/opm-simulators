/*
  Copyright 2021 Equinor.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/simulators/utils/DamarisKeywords.hpp>
#include <string>
#include <map>

/*
    Below is the Damaris Keywords supported by Damaris to be filled
    in the built-in XML file.

    The entries in the map below will be filled by the corresponding
    Damaris Keywords. Yet, only output directory and FileMode are to
    be chosen by the user
*/

namespace Opm::DamarisOutput
{
std::map<std::string, std::string>
DamarisKeywords(std::string OutputDir, 
                    bool enableDamarisOutputCollective, 
                    bool saveToHDF5, 
                    int  nDamarisCores,
                    int  nDamarisNodes,
                    long shmemSizeBytes,
                    std::string pythonFilename, 
                    std::string simName, 
                    std::string logLevel,
                    std::string paraviewPythonFilename )
{
    std::string saveToHDF5_str("MyStore") ;
    if (! saveToHDF5 )  saveToHDF5 = "" ;
    
    std::string publishToPython_str("") ;
    if (pythonFilename != "") publishToPython_str="PythonScript" ; // the name of the PyScript XML element
    
    std::string damarisOutputCollective_str("") ;
    if (enableDamarisOutputCollective) {
        damarisOutputCollective_str="Collective"
    } else {
        damarisOutputCollective_str="FilePerCore"
    }
    
    std::string  simName_str  ;
    if (simName == "") {
        // Having a different simulation name is important if multiple simulations 
        // are running on the same node, as it is used to name the simulations shmem area
        // and when one sim finishes it removes its shmem file.
        simName_str = "opm-sim-" + damaris::Environment::GetMagicNumber() ;
    } else {
        simName_str = simName ;
    }
    
    if (nDamarisCores > 0) && (nDamarisNodes > 0)
    {
        nDamarisNodes = 0 ; // Default is to use Damaris Cores
    }
    string nDamarisCores_str  ;
    if ( nDamarisCores != 0 ) {
        nDamarisCores_str = to_string(nDamarisCores);
    } else {
        nDamarisCores_str = "0" ;
    }
    
    string nDamarisNodes_str  ;
    if ( nDamarisNodes_str != 0 ) {
        nDamarisNodes_str = to_string(nDamarisNodes);
    } else {
        nDamarisNodes_str = "0" ;
    }
    
    string shmemSizeBytes_str  ;
    if ( shmemSizeBytes != 0 ) {
        shmemSizeBytes_str = to_string(shmemSizeBytes);
    } else {
        shmemSizeBytes_str = "536870912" ;
    }

   // _MAKE_AVAILABLE_IN_PYTHON_
   // _PYTHON_SCRIPT_

    std::map<std::string, std::string> damaris_keywords = {
        {"_SHMEM_BUFFER_BYTES_REGEX_", shmemSizeBytes_str},
        {"_DC_REGEX_", nDamarisCores_str},
        {"_DN_REGEX_", nDamarisNodes_str},
        {"_File_Mode", damarisOutputCollective_str},
        {"_MORE_VARIABLES_REGEX_", ""},
        {"_PATH_REGEX_", OutputDir},
        {"_MYSTORE_OR_EMPTY_REGEX_", saveToHDF5_str},
        {"_PARAVIEW_PYTHON_SCRIPT_",""},  /* this has to be before _PYTHON_SCRIPT_ entry */
        {"_PYTHON_SCRIPT_",pythonFilename}, /* if a Python script is speified then assume that we want to publish the data to Python */
        {"_PRESSURE_UNIT_","Pa"},
        {"_MAKE_AVAILABLE_IN_PYTHON_",publishToPython_str},  /* must match  <pyscript name="PythonScript" */
        {"_SIM_NAME_",simName_str},
        {"_LOG_LEVEL_",logLevel_str},
    };
    return damaris_keywords;
    
}

} // namespace Opm::DamarisOutput
