#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import datetime

import pcraster as pcr
from pcraster.framework import DynamicModel
from pcraster.framework import DynamicFramework

# time object
from currTimeStep import ModelTime

from outputNetcdf import OutputNetcdf
import virtualOS as vos

import logging
logger = logging.getLogger(__name__)

class CalcFramework(DynamicModel):

    def __init__(self, cloneMapFileName,\
                       modelTime, \
                       input_files, \
                       output_files
                       ):
        DynamicModel.__init__(self)
        
        # set the clone map
        self.cloneMapFileName = cloneMapFileName
        pcr.setclone(self.cloneMapFileName)
        
        # time variable/object
        self.modelTime = modelTime
        
        # a dictionary containing input files
        self.input_files = input_files
        
        # a dictionary containing output files
        self.output_files  = output_files 
        self.output_folder = self.output_files["folder"]
        
        # prepare temporary directory
        self.tmpDir = self.output_folder + "/tmp/"
        try:
            os.makedirs(self.tmpDir)
        except:
            os.system('rm -r '+tmpDir+"/*")
        
        # cell area (m2)
        self.cell_area_total = vos.readPCRmapClone(v = self.input_files["cell_area"], \
                                                   cloneMapFileName = self.cloneMapFileName, \
                                                   tmpDir = self.tmpDir
                                                   )
        
        # ldd                                            
        ldd = vos.readPCRmapClone(v = self.input_files["ldd_map"], \
                                  cloneMapFileName = self.cloneMapFileName, \
                                  tmpDir = self.tmpDir
                                  absolutePath = None, isLddMap = True, cover = None, isNomMap = False
                                  )
        self.ldd = pcr.lddrepair(pcr.ldd(ldd))                                   
        
        # object for reporting
        self.netcdf_report = OutputNetcdf(mapattr_dict = None,\
                                          cloneMapFileName = cloneMapFileName,\
                                          netcdf_format = "NETCDF4",\
                                          netcdf_zlib = False,\
                                          netcdf_attribute_dict = None,)       

        
    def initial(self): 

        # general attributes for netcdf output files
        attributeDictionary = {}
        attributeDictionary['institution']   = "Department of Physical Geography, Utrecht University, the Netherlands"
        attributeDictionary['history'    ]   = "Files are created on " + str(datetime.datetime.now())
        attributeDictionary['references' ]   = "See description."
        attributeDictionary['source'     ]   = "See description."
        attributeDictionary['comment'    ]   = "See description. Calculated on the folder " + str(self.output_files["folder"]) 
        attributeDictionary['disclaimer' ]   = "Great care was exerted to prepare these data. Notwithstanding, use of the model and/or its outcome is the sole responsibility of the user." 

        # make a netcdf output file for monthly estimate irrigation demand
        attributeDictionary['title'      ]   = "XXX."
        attributeDictionary['description']   = "XXX. "
        self.netcdf_report.createNetCDF(self.output_files["gw_discharge_contribution"],\
                                        "gw_discharge_contribution",\
                                        "-",\
                                        "gw_discharge_contribution",\
                                        attributeDictionary)

    def dynamic(self):
        
        # re-calculate current model time using current pcraster timestep value
        self.modelTime.update(self.currentTimeStep())
        
        # calculate gw_discharge_contribution
        if self.modelTime.isLastDayOfMonth():

            # read input files
            # - groundwater discharge (m3/day)
            local_gw_discharge = (1/self.modelTime.day) * self.cell_area_total *\
                                 vos.netcdf2PCRobjClone(ncFile            = self.input_files["local_gw_discharge"],\
                                                        varName           = "automatic",\
                                                        dateInput         = self.modelTime.fulldate,\
                                                        useDoy            = None,\
                                                        cloneMapFileName  = self.cloneMapFileName)
            # - interflow (m3/day)
            interflow          = (1/self.modelTime.day) * self.cell_area_total *\
                                 vos.netcdf2PCRobjClone(ncFile            = self.input_files["interflow"],\
                                                        varName           = "automatic",\
                                                        dateInput         = self.modelTime.fulldate,\
                                                        useDoy            = None,\
                                                        cloneMapFileName  = self.cloneMapFileName)
            # - direct runoff (m3/day)
            direct_runoff      = (1/self.modelTime.day) * self.cell_area_total *\
                                 vos.netcdf2PCRobjClone(ncFile            = self.input_files["direct_runoff"],\
                                                        varName           = "automatic",\
                                                        dateInput         = self.modelTime.fulldate,\
                                                        useDoy            = None,\
                                                        cloneMapFileName  = self.cloneMapFileName)
            
            
            # shall we consider negative values of discharge?
            local_gw_discharge = pcr.max(0.0, local_gw_discharge)
            
            # accumulating direct runoff, interflow and groundwater discharge/baseflow - unit: m3/day
            accumulated_runoff = pcr.catchmenttotal( (direct_runoff + interflow + local_gw_discharge), self.ldd)
            
            # percentage contribution of gw_discharge
            gw_discharge_contribution = local_gw_discharge / (accumulated_runoff - direct_runoff + interflow + local_gw_discharge)
            

        # save monthly irrigation demand and monthly irrigation requirement to files (km3/month)
        if self.modelTime.isLastDayOfMonth():

            # reporting
            timeStamp = datetime.datetime(self.modelTime.year,\
                                          self.modelTime.month,\
                                          self.modelTime.day,0)
            varFields = {}
            varFields["gw_discharge_contribution"] = pcr.pcr2numpy(gw_discharge_contribution, vos.MV)
            self.netcdf_report.dataList2NetCDF(self.output_files["estimate_irrigation_demand"],\
                                               ["gw_discharge_contribution"],\
                                               varFields,\
                                               timeStamp)


def main():
    
    # ~ # use the following system arguments
    # ~ start_year                          = sys.argv[1]
    # ~ end_year                            = sys.argv[2]
    # ~ pcrglobwb_input_folder              = sys.argv[3]
    # ~ irrigated_area_in_hectar_input_file = sys.argv[4]
    # ~ pcrglobwb_monthly_output_folder     = sys.argv[5]
    # ~ pcrglobwb_daily_output_folder       = sys.argv[6]
    # ~ output_folder_for_irrigation_demand = sys.argv[7]
    # ~ output_file_for_irrigation_demand   = sys.argv[8]
    
    # ~ # starting and end date
    # ~ startDate = "%s-01-01" % (str(start_year))
    # ~ endDate   = "%s-12-31" % (str(end_year))

    startDate = "1979-01-01"
    endDate   = "2013-12-31"
    # - TODO: Using the endDate = 2019-12-31 (after the MODFLOW run is done).

    # a dictionary containing input files
    input_files = {}
    
    # input from PCR-GLOBWB and MODFLOW INPUT files (model parameters) - at 30sec resolution
    # - clone map (-), cell area (m2)
    input_files["clone_map"]                    = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/cloneMaps/australia_30sec.map"    
    input_files["ldd_map"]                      = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/routing/surface_water_bodies/version_2020-05-XX/lddsound_30sec_version_202005XX_correct_lat.nc"
    input_files["cell_area"]                    = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/others/estimate_cell_dimension/30sec/cdo_grid_area_30sec_map_correct_lat.nc"


    # input from PCR-GLOBWB run OUTPUT files - at 5 arcmin resolution
    # - direct runoff and interflow (m.month-1) 
    input_files["direct_runoff"]      = "/scratch/depfg/otoo0001/test_australia_w5e5/netcdf/directRunoff_monthTot_output.nc"
    input_files["interflow"]          = "/scratch/depfg/otoo0001/test_australia_w5e5/netcdf/interflowTotal_monthTot_output.nc"
    # - TODO: Perhaps, we also should think other components, such as surface water abstraction and evaporation from surface water


    # input from MODFLOW run OUTPUT files - at 30 sec resolution
    input_files["local_gw_discharge"] = "/scratch/depfg/sutan101/test_pcrglobwb_gmglob_30sec/australia_test_satareafrac_w5e5/australia_30sec/transient/netcdf/baseflow_monthEnd_output.nc" 


    # a dictionary containing output files
    output_files = {}
    output_files["folder"]                    = "/scratch/depfg/sutan101/test_gw_discharge_contribution/"
    output_files["gw_discharge_contribution"] = output_files["folder"] + "/local_gw_discharge_contribution.nc" 


    # make output folder
    output_folder = output_files["folder"]
    try:
        os.makedirs(output_folder)
    except:
        os.system('rm -r ' + output_folder + "/*") ### THIS IS DANGEROUS
        pass


    # prepare logger and its directory
    log_file_location = output_folder + "/log/"
    try:
        os.makedirs(log_file_location)
    except:
        cmd = 'rm -r ' + log_file_location + "/*"
        os.system(cmd)
        pass
    vos.initialize_logging(log_file_location)
    

    # time object
    modelTime = ModelTime() # timeStep info: year, month, day, doy, hour, etc
    modelTime.getStartEndTimeSteps(startDate, endDate)
    
    
    # put all of them in the pcraster dynamic calculation framework
    calculationModel = CalcFramework(input_files["clone_map"],\
                                     modelTime, \
                                     input_files, \
                                     output_files)

    dynamic_framework = DynamicFramework(calculationModel, modelTime.nrOfTimeSteps)
    dynamic_framework.setQuiet(True)
    dynamic_framework.run()


if __name__ == '__main__':
    sys.exit(main())
