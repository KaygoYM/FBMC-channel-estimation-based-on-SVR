%author:KAI
function [BinaryDataStream_FBMC_Aux,x_FBMC_Aux,s_FBMC_Aux]= FBMC_data(AuxiliaryMethod,FBMC,PAM)
    BinaryDataStream_FBMC_Aux = randi([0 1],AuxiliaryMethod.NrTransmittedSymbols*log2(PAM.ModulationOrder),1);
    xD_FBMC_Aux = PAM.Bit2Symbol(BinaryDataStream_FBMC_Aux);
    x_FBMC_Aux = reshape(xD_FBMC_Aux,FBMC.Nr.Subcarriers,FBMC.Nr.MCSymbols);
    s_FBMC_Aux = FBMC.Modulation(x_FBMC_Aux);
end
