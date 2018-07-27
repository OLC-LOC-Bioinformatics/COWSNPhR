import os


class Get_Specie_Parameters_Step1():

    def __init__(self):
        real_path = os.path.dirname(os.path.realpath(__file__))
        print("real path command --> {}".format(real_path))
        real_path = real_path.split('/')
        root_path = '/'.join(real_path)
        self.dependents_dir = root_path + "/dependencies"
        if os.path.isdir("/bioinfo11/TStuber/Results"): #check bioinfo from server
            self.upload_to = "/bioinfo11/TStuber/Results"
        else:
            self.upload_to = None

    def choose(self, species_selection):
        if species_selection == "salmonella":
            script_dependents = str(self.dependents_dir) + "/gen-bact/salmonella/snp_pipeline/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to),
                "spoligo_db": None,
                "reference": script_dependents + "/NC_016856-NC_016855.fasta",
                "hqs": script_dependents + "/NC_016856-NC_016855HighestQualitySNPs.vcf",
                "gbk_file": script_dependents + "/NC_016856-NC_016855.gbk",
                "species": species_selection,
            }
            return(parameters)

        elif species_selection == "ab1":
            script_dependents = str(self.dependents_dir) + "/brucella/abortus1/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/abortus1/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_00693c.fasta",
                "hqs": script_dependents + "/NC_00693cHighestQualitySNPs.vcf",
                "gbk_file": script_dependents + "/NC_006932-NC_006933.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "ab3":
            script_dependents = str(self.dependents_dir) + "/brucella/abortus3/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/abortus3/data",
                "spoligo_db": None,
                "reference": script_dependents + "/CP007682-7683c.fasta",
                "hqs": script_dependents + "/CP007682-7683cHighestQualitySNPs.vcf",
                "gbk_file": script_dependents + "/CP007682-CP007683.gbk",
                "species": species_selection,
            }
            return (parameters)
        elif species_selection == "canis":
            script_dependents = str(self.dependents_dir) + "/brucella/canis/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/canis/data",
                "spoligo_db": None,
                "reference": script_dependents + "/BcanisATCC23365.fasta",
                "hqs": script_dependents + "/canisHighestQualitySNPs.vcf",
                "gbk_file": script_dependents + "/NC_010103-NC_010104.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "ceti1":
            script_dependents = str(self.dependents_dir) + "/brucella/ceti1/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/ceti1/data",
                "spoligo_db": None,
                "reference": script_dependents + "/Bceti1Cudo.fasta",
                "hqs": script_dependents + "/ceti1HighestQualitySNPs.vcf",
                "gbk_file": None,
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "ceti2":
            script_dependents = str(self.dependents_dir) + "/brucella/ceti2/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/ceti2/data",
                "spoligo_db": None,
                "reference": script_dependents + "/Bceti2-TE10759.fasta",
                "hqs": script_dependents + "/ceti2HighestQualitySNPs.vcf",
                "gbk_file": script_dependents + "/NC_022905-NC_022906.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "mel1":
            script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv1/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/melitensis-bv1/data",
                "spoligo_db": None,
                "reference": script_dependents + "/mel-bv1-NC003317.fasta",
                "hqs": script_dependents + "/mel-bv1-NC003317-highqualitysnps.vcf",
                "gbk_file": script_dependents + "/NC_003317-NC_003318.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "mel1b":
            script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv1b/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/melitensis-bv1b/data",
                "spoligo_db": None,
                "reference": script_dependents + "/mel-bv1b-CP018508.fasta",
                "hqs": script_dependents + "/mel-bv1b-CP018508-highqualitysnps.vcf",
                "gbk_file": script_dependents + "/mel-bv1b-CP018508.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "mel2":
            script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv2/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/melitensis-bv2/data",
                "spoligo_db": None,
                "reference": script_dependents + "/mel-bv2-NC012441.fasta",
                "hqs": script_dependents + "/mel-bv2-NC012441-highqualitysnps.vcf",
                "gbk_file": script_dependents + "/NC_012441-NC_012442.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "mel3":
            script_dependents = str(self.dependents_dir) + "/brucella/melitensis-bv3/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/melitensis-bv3/data",
                "spoligo_db": None,
                "reference": script_dependents + "/mel-bv3-NZCP007760.fasta",
                "hqs": script_dependents + "/mel-bv3-NZCP007760-highqualitysnps.vcf",
                "gbk_file": script_dependents + "/NZ_CP007760-NZ_CP007761.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "suis1":
            script_dependents = str(self.dependents_dir) + "/brucella/suis1/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/suis1/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_017251-NC_017250.fasta",
                "hqs": script_dependents + "/B00-0468-highqualitysnps.vcf",
                "gbk_file": script_dependents + "/NC_017251-NC_017250.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "suis2":
            script_dependents = str(self.dependents_dir) + "/brucella/suis2/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/suis2/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_010169-NC_010167.fasta",
                "hqs": script_dependents + "/BsuisF7-06-1-highqualitysnps.vcf",
                "gbk_file": script_dependents + "/NC_010169-NC_010167.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "suis3":
            script_dependents = str(self.dependents_dir) + "/brucella/suis3/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/suis3/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NZ_CP007719-NZ_CP007718.fasta",
                "hqs": script_dependents + "/highqualitysnps.vcf",
                "gbk_file": None,
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "suis4":
            script_dependents = str(self.dependents_dir) + "/brucella/suis4/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/suis4/data",
                "spoligo_db": None,
                "reference": script_dependents + "/B-REF-BS4-40.fasta",
                "hqs": script_dependents + "/suis4HighestQualitySNPs.vcf",
                "gbk_file": None,
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "ovis":
            script_dependents = str(self.dependents_dir) + "/brucella/ovis/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/ovis/data",
                "spoligo_db": None,
                "reference": script_dependents + "/BovisATCC25840.fasta",
                "hqs": script_dependents + "/BovisATCC25840HighestQualitySNPs.vcf",
                "gbk_file": script_dependents + "/NC_009505-NC_009504.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "neo":
            script_dependents = str(self.dependents_dir) + "/brucella/neotomae/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/brucella/neotomae/data",
                "spoligo_db": None,
                "reference": script_dependents + "/KN046827.fasta",
                "hqs": script_dependents + "/ERR1845155-highqualitysnps.vcf",
                "gbk_file": script_dependents + "/KN046827.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "af":
            script_dependents = str(self.dependents_dir) + "/mycobacterium/tbc/af2122/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacterium/tbc/af2122/script1",
                "spoligo_db": script_dependents + "/spoligotype_db.txt",
                "reference": script_dependents + "/NC_002945v4.fasta",
                "hqs": script_dependents + "/highqualitysnps.vcf",
                "gbk_file": script_dependents + "/NC_002945v4.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "h37":
            script_dependents = str(self.dependents_dir) + "/mycobacterium/tbc/h37/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacterium/tbc/h37/script1",
                "spoligo_db": script_dependents + "/spoligotype_db.txt",
                "reference": script_dependents + "/NC000962.fasta",
                "hqs": script_dependents + "/15-3162-highqualitysnps.vcf",
                "gbk_file": script_dependents + "/NC_000962.gbk",
                "species": species_selection,
            }
            return(parameters)
        elif species_selection == "para":
            script_dependents = str(self.dependents_dir) + "/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script1"
            parameters = {
                "upload_to": str(self.upload_to) + "/mycobacterium/avium_complex/para_cattle-bison/data",
                "spoligo_db": None,
                "reference": script_dependents + "/NC_002944.fasta",
                "hqs": script_dependents + "/HQ-NC002944.vcf",
                "gbk_file": script_dependents + "/NC_002944.gbk",
                "species": species_selection,
            }
            return (parameters)
        else:
            parameters = {
                "upload_to": None,
                "spoligo_db": None,
                "reference": None,
                "hqs": None,
                "gbk_file": None,
                "species": None,
            }
            return (parameters)


class Get_Specie_Parameters_Step2():

    def __init__(self):
        real_path = os.path.dirname(os.path.realpath(__file__))
        print("real path command --> {}".format(real_path))
        real_path = real_path.split('/')
        root_path = '/'.join(real_path)
        self.dependents_dir = root_path + "/dependencies"
        if os.path.isdir("/bioinfo11/TStuber/Results"): # check bioinfo from server
            self.upload_to = "/bioinfo11/TStuber/Results"
        else:
            self.upload_to = None
    ###### bruc_private_codes(upload_to)
    def choose(self, species_selection):
        if species_selection == "af":
            script_dependents = self.dependents_dir + "/mycobacterium/tbc/af2122/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 150,
                "N_gatk_threshold": 150,
                "genotypingcodes": str(self.upload_to) + "/mycobacterium/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_002945v4.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/Filtered_Regions.xlsx", # previous excelinfile
                "step2_upload": str(self.upload_to) + "/mycobacterium/tbc/af2122/script2", #previous bioinfoVCF
            }
        elif species_selection == "salmonella":
            script_dependents = self.dependents_dir + "/gen-bact/salmonella/snp_pipeline/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": None,
                "gbk_file": script_dependents + "/NC_016856-NC_016855.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx"
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx"
                "filter_file": script_dependents + "/filter_files"
                "step2_upload": str(self.upload_to) + "/gen-bact/salmonella/snp_pipeline/script2",
            }
        elif species_selection == "suis1":
            script_dependents = self.dependents_dir + "/brucella/suis1/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_017251-NC_017250.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/suis1/vcfs",
            }
        elif species_selection == "suis2":
            script_dependents = self.dependents_dir + "/brucella/suis2/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_010169-NC_010167.gbk",
                "definingSNPs": script_dependents + "/Defining_SNPs.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/suis2/vcfs",
            }
        elif species_selection == "suis3":
            script_dependents = self.dependents_dir + "/brucella/suis3/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NZ_CP007719-NZ_CP007718.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/suis3/vcfs",
            }
        elif species_selection == "suis4":
            script_dependents = self.dependents_dir + "/brucella/suis4/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": None,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files"
                "step2_upload": str(self.upload_to) + "/brucella/suis4/vcfs",
            }
        elif species_selection == "ab1":
            script_dependents = self.dependents_dir + "/brucella/abortus1/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_006932-NC_006933.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/abortus1/vcfs",
            }
        elif species_selection == "ab3":
            script_dependents = self.dependents_dir + "/brucella/abortus3/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/CP007682-CP007683.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/abortus3/vcfs",
            }
        elif species_selection == "mel1":
            script_dependents = self.dependents_dir + "/brucella/melitensis-bv1/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_003317-NC_003318.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/melitensis-bv1/vcfs",
            }
        elif species_selection == "mel1b":
            script_dependents = self.dependents_dir + "/brucella/melitensis-bv1b/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/mel-bv1b-CP018508.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/melitensis-bv1b/vcfs",
            }
        elif species_selection == "mel2":
            script_dependents = self.dependents_dir + "/brucella/melitensis-bv2/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_012441-NC_012442.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/melitensis-bv2/vcfs",
            }
        elif species_selection == "mel3":
            script_dependents = self.dependents_dir + "/brucella/melitensis-bv3/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NZ_CP007760-NZ_CP007761.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/melitensis-bv3/vcfs",
            }
        elif species_selection == "canis":
            script_dependents = self.dependents_dir + "/brucella/canis/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_010103-NC_010104.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/canis/vcfs",
            }
        elif species_selection == "ceti1":
            script_dependents = self.dependents_dir + "/brucella/ceti1/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": None,
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/ceti1/vcfs",
            }
        elif species_selection == "ceti2":
            script_dependents = self.dependents_dir + "/brucella/ceti2/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_022905-NC_022906.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/ceti2/vcfs",
            }   
        elif species_selection == "ovis":
            script_dependents = self.dependents_dir + "/brucella/ovis/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_009505-NC_009504.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/ovis/vcfs",
            }     
        elif species_selection == "neo":
            script_dependents = self.dependents_dir + "/brucella/neotomae/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 300,
                "N_gatk_threshold": 350,
                "genotypingcodes": str(self.upload_to) + "/brucella/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/KN046827.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/brucella/neotomae/vcfs",
            }
        elif species_selection == "bovis":
            script_dependents = self.dependents_dir + "/mycobacterium/tbc/tbbov/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 150,
                "N_gatk_threshold": 150  ,
                "genotypingcodes": str(self.upload_to) + "/mycobacterium/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_002945.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/mycobacterium/tbc/tbbov/script2",
            }
        elif species_selection == "af":
            script_dependents = self.dependents_dir + "/mycobacterium/tbc/af2122/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 150,
                "N_gatk_threshold": 150 ,
                "genotypingcodes": str(self.upload_to) + "/mycobacterium/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_002945v4.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/mycobacterium/tbc/af2122/script2",
            }
        elif species_selection == "h37":
            script_dependents = self.dependents_dir + "/mycobacterium/tbc/h37/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 150,
                "N_gatk_threshold": 150,
                "genotypingcodes": str(self.upload_to) + "/mycobacterium/genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_000962.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations_python.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/mycobacterium/tbc/h37/script2",
            }
        elif species_selection == "para":
            script_dependents = self.dependents_dir + "/mycobacterium/avium_complex/para_cattle-bison/script_dependents/script2"
            parameters = {
                "qual_gatk_threshold": 150,
                "N_gatk_threshold": 150,
                "genotypingcodes": str(self.upload_to) + "/mycobacterium/avium_genotyping_codes.xlsx",
                "gbk_file": script_dependents + "/NC_002944.gbk",
                "definingSNPs": script_dependents + "/DefiningSNPsGroupDesignations.xlsx",
                "remove_from_analysis": script_dependents + "/RemoveFromAnalysis.xlsx",
                "filter_file": script_dependents + "/filter_files",
                "step2_upload": str(self.upload_to) + "/mycobacterium/avium_complex/para_cattle-bison/vcfs",
            }
        else:
            parser.print_help()
            print ("\n#####EXIT AT SETTING OPTIONS, Check that a \"-s\" species was provided\n")
            sys.exit(0)

        if self.upload_to is None:
            parameters["genotypingcodes"] = None
            parameters["step2_upload"] = None
        return(parameters)