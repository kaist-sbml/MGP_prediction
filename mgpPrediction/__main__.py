from mgpPrediction import inputParser
from mgpPrediction.statMGPprediction import runAnalysis
from mgpPrediction.standardizedComp import standardizedCompare
from mgpPrediction.utils import createOrganizedOutput


def main():
    parser = inputParser.argument_parser()
    options = parser.parse_args()
    
    output_dir = options.output_dir
    flux_path = options.flux
    mutation_path = options.mutation
    
    averagedTargetFluxsum, pValueINFO=runAnalysis(flux_path, mutation_path)
    MGPs = standardizedCompare(averagedTargetFluxsum, pValueINFO)
    
    createOrganizedOutput(MGPs).to_csv(output_dir+'/predictedMGPs.csv', index=False)