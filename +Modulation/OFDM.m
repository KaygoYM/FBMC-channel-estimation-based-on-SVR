classdef OFDM < handle 
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2016 by Institute of Telecommunications, TU Wien
% www.tc.tuwien.ac.at    


   properties (SetAccess = private)
       Nr
       PHY
       Implementation
   end
   
   
   methods
      %% Class constructor, define default values. 
      function obj = OFDM(varargin)
         %% Initialize parameters, set default values
         if numel(varargin)==8
         	obj.Nr.Subcarriers = varargin{1};
         	obj.Nr.MCSymbols = varargin{2};  
            obj.PHY.SubcarrierSpacing = varargin{3};
            obj.PHY.SamplingRate = varargin{4};
            obj.PHY.IntermediateFrequency = varargin{5};
            obj.PHY.TransmitRealSignal = varargin{6};
            obj.PHY.CyclicPrefixLength = varargin{7};
            obj.PHY.ZeroGuardTimeLength = varargin{8};
         elseif numel(varargin)==0
            % LTE Default Values (can later be changed using FBMC.Set...)
            obj.Nr.Subcarriers = 12;
            obj.Nr.MCSymbols = 14; 
            obj.PHY.SubcarrierSpacing = 15e3; 
            obj.PHY.SamplingRate = obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing; 
            obj.PHY.IntermediateFrequency = 0; 
            obj.PHY.TransmitRealSignal = false; 
            obj.PHY.CyclicPrefixLength = ceil(4.7e-6*obj.PHY.SamplingRate)/obj.PHY.SamplingRate;
            obj.PHY.ZeroGuardTimeLength = 0;
         else 
            error('Number of input variables must be either 0 (default values) or 10');
         end
         
         %% Check Parameters
         if mod(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing),1)~=0
            obj.PHY.SubcarrierSpacing=obj.PHY.SamplingRate/(round(obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)));         
            disp('Sampling Rate must be a multiple of the subcarrier spacing!');
            disp(['Therefore, the subcarrier spacing is set to: ' int2str(obj.PHY.SubcarrierSpacing) 'Hz']);            
         end
         
         if mod(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing,2)~=0
            obj.PHY.IntermediateFrequency = round(obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing)*obj.PHY.SubcarrierSpacing;
            disp('The intermediate frequency must be a multiple of the subcarrier spacing!');
            disp(['Therefore, the intermediate frequency is set to ' int2str(obj.PHY.IntermediateFrequency) 'Hz']);            
         end
         
         if (obj.PHY.SamplingRate<obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing)
             error('Sampling Rate must be higher: at least Number of Subcarriers times Subcarrier Spacing');
         end
         
         if abs(mod(obj.PHY.CyclicPrefixLength*obj.PHY.SamplingRate,1))>10^(-14)
             obj.PHY.CyclicPrefixLength=round(obj.PHY.CyclicPrefixLength*obj.PHY.SamplingRate)/obj.PHY.SamplingRate;
             disp('The length of the cyclic prefix times the sampling rate must be a integer!');
             disp(['Therefore, the cyclic prefix length is set to ' num2str(obj.PHY.CyclicPrefixLength) 's']);         
         end
         obj.Implementation.CyclicPrefix = round(obj.PHY.CyclicPrefixLength*obj.PHY.SamplingRate);
         obj.Implementation.ZeroGuardSamples = round(obj.PHY.ZeroGuardTimeLength*obj.PHY.SamplingRate);
         
         %% Dependent parameters
         obj.PHY.dt = 1/obj.PHY.SamplingRate;
            
         obj.Implementation.TimeSpacing = obj.PHY.SamplingRate/(obj.PHY.SubcarrierSpacing)+obj.Implementation.CyclicPrefix; 
         obj.PHY.TimeSpacing = obj.Implementation.TimeSpacing*obj.PHY.dt;
          
         obj.Nr.SamplesTotal = (obj.Nr.MCSymbols*obj.Implementation.TimeSpacing)+2*obj.Implementation.ZeroGuardSamples;    
         
         obj.Implementation.FFTSize = obj.PHY.SamplingRate/obj.PHY.SubcarrierSpacing;

         obj.Implementation.IntermediateFrequency = obj.PHY.IntermediateFrequency/obj.PHY.SubcarrierSpacing;
         
         % Normalization factor so that power = 1       
         obj.Implementation.NormalizationFactor = sqrt(obj.PHY.SamplingRate^2/obj.PHY.SubcarrierSpacing^2/obj.Nr.Subcarriers);    
      
      end
      
      %% Set Functions
      function SetNrSubcarriers(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.Nr.Subcarriers = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end      
      function SetNrMCSymbols(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.Nr.MCSymbols = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end
      function SetSubcarrierSpacing(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.PHY.SubcarrierSpacing = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end 
      function SetSamplingRate(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.PHY.SamplingRate = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end       
      function SetIntermediateFrequency(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.PHY.IntermediateFrequency = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end        
      function SetTransmitRealSignal(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.PHY.TransmitRealSignal = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end         
      function SetMethod(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.Method = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end              
      function SetOverlappingFactor(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.PrototypeFilter.OverlappingFactor = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end                  
      function SetInitialPhaseShift(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.Implementation.InitialPhaseShift = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end                  
      function SetUsePolyphase(varargin) % This is a very sloppy implementation  :/
         obj = varargin{1}; obj.Implementation.UsePolyphase = varargin{2};
         NEWobj = Modulation.FBMC(obj.Nr.Subcarriers,obj.Nr.MCSymbols,obj.PHY.SubcarrierSpacing,obj.PHY.SamplingRate,obj.PHY.IntermediateFrequency,obj.PHY.TransmitRealSignal,obj.Method,obj.PrototypeFilter.OverlappingFactor,obj.Implementation.InitialPhaseShift,obj.Implementation.UsePolyphase);
         obj.Method = NEWobj.Method; obj.Nr =  NEWobj.Nr; obj.PHY =  NEWobj.PHY;  obj.PrototypeFilter =  NEWobj.PrototypeFilter;  obj.Implementation =  NEWobj.Implementation;
      end                
                
      %% Modulation and Demodulation
      function TransmitSignal = Modulation(varargin)
          %asdf
         obj = varargin{1};
         DataSymbols = varargin{2};
                  
         DataSymbolsTemp = zeros(obj.Implementation.FFTSize,obj.Nr.MCSymbols); 
         DataSymbolsTemp(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers),:) = DataSymbols.*obj.Implementation.NormalizationFactor;
         if obj.PHY.TransmitRealSignal
             DataSymbolsTemp = (DataSymbolsTemp+conj(DataSymbolsTemp([1 end:-1:2],:)))/sqrt(2);
         end 
         TransmitSignalNoCP =  ifft(DataSymbolsTemp);
         TransmitSignal =[zeros(obj.Implementation.ZeroGuardSamples,1);reshape([TransmitSignalNoCP(end-obj.Implementation.CyclicPrefix+1:end,:);TransmitSignalNoCP],obj.Implementation.TimeSpacing*obj.Nr.MCSymbols,1);zeros(obj.Implementation.ZeroGuardSamples,1)];               
      end     
      function ReceivedSymbols = Demodulation(varargin)
         obj = varargin{1};
         ReceivedSignal = varargin{2};
         
         ReceivedSignal_reshape = reshape(ReceivedSignal((obj.Implementation.ZeroGuardSamples+1):(end-obj.Implementation.ZeroGuardSamples)),obj.Implementation.TimeSpacing,obj.Nr.MCSymbols);
         ReceivedSymbolsTemp = fft(ReceivedSignal_reshape(obj.Implementation.CyclicPrefix+1:end,:));
         if obj.PHY.TransmitRealSignal
%             ReceivedSymbolsTemp = (ReceivedSymbolsTemp+conj(ReceivedSymbolsTemp([1 end:-1:2],:)))/sqrt(2);
            ReceivedSymbolsTemp = ReceivedSymbolsTemp*sqrt(2);
         end 
         
         ReceivedSymbols = ReceivedSymbolsTemp(obj.Implementation.IntermediateFrequency+(1:obj.Nr.Subcarriers),:)/(obj.Implementation.NormalizationFactor);      
      end
            
      %% Matrix Description
      function TXMatrix = GetTXMatrix(obj) 
         if obj.PHY.TransmitRealSignal
            % Does nor work because conjungate complex is not a linear operation for complex symols. For FBMC it worked because real symbols are transmitted.
            error('GetTXMatrix is not supported for PHY.TransmitRealSignal == true!')        
         end    
         TXMatrix=zeros(obj.Nr.SamplesTotal,obj.Nr.Subcarriers*obj.Nr.MCSymbols);
         TXMatrixTemp=zeros(obj.Nr.SamplesTotal,obj.Nr.Subcarriers);
         x = zeros(obj.Nr.Subcarriers, obj.Nr.MCSymbols);
         for i_l= 1:obj.Nr.Subcarriers;
            x(i_l)=1;
            TXMatrixTemp(:,i_l) = obj.Modulation(x);
            x(i_l)=0;
         end
         for i_k=1:obj.Nr.MCSymbols
            TXMatrix(:,(1:obj.Nr.Subcarriers)+(i_k-1)*obj.Nr.Subcarriers)=circshift(TXMatrixTemp,[(i_k-1)*obj.Implementation.TimeSpacing,0]);           
         end                       
      end
      function RXMatrix = GetRXMatrix(obj) 
         if obj.PHY.TransmitRealSignal
            obj.PHY.TransmitRealSignal = false;
            RXMatrix = sqrt(2)*obj.GetTXMatrix'*(obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing/(obj.PHY.SamplingRate));             
            obj.PHY.TransmitRealSignal = true;
         else
            RXMatrix = obj.GetTXMatrix'*(obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing/(obj.PHY.SamplingRate));            
         end
         IndexTemp = obj.Implementation.ZeroGuardSamples+bsxfun(@plus, (1:obj.Implementation.CyclicPrefix)', (0:obj.Nr.MCSymbols-1)*obj.Implementation.TimeSpacing);
         RXMatrix(:,IndexTemp(:))=0;
      end
                              
      %% Plot
      function [TransmitPower,Time]=PlotTransmitPower(varargin)
          obj = varargin{1};
          if numel(varargin)==2
            [V,D] = eig(varargin{2});
          else
            V = eye(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
            D = V;
          end         
          D=sqrt(D);
          TransmitPower = zeros(obj.Nr.SamplesTotal,1);
          for i_lk = 1:obj.Nr.Subcarriers*obj.Nr.MCSymbols    
%             TransmitPower = TransmitPower+abs(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2;        
            % The modulation is not linear => we have to consider 1 and +j.
            % Matters only for OFDM.PHY.TransmitRealValuedSignal == 0 ;
            % Still should be checked when time!
            TransmitPower = TransmitPower+(abs(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2+abs(obj.Modulation(reshape(j*V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk)).^2)/2;
            if mod(i_lk,1000)==0
                disp([int2str(i_lk/(obj.Nr.Subcarriers*obj.Nr.MCSymbols )*100) '%']);
            end 
          end
          Time = (0:length(TransmitPower)-1)*obj.PHY.dt;
          if nargout==0
            plot(Time,TransmitPower);
            ylabel('Transmit Power');
            xlabel('Time(s)'); 
          end
      end        
      function [PowerSpectralDensity,Frequency] = PlotPowerSpectralDensity(varargin)
        obj = varargin{1};
        if numel(varargin)==2
            [V,D] = eig(varargin{2});
        else
            V = eye(obj.Nr.Subcarriers*obj.Nr.MCSymbols);
            D = V;
        end         
        D=sqrt(D);
        PowerSpectralDensity = zeros(obj.Nr.SamplesTotal,1);
        for i_lk = 1:obj.Nr.Subcarriers*obj.Nr.MCSymbols    
            PowerSpectralDensity = PowerSpectralDensity+abs(fft(obj.Modulation(reshape(V(:,i_lk),obj.Nr.Subcarriers,obj.Nr.MCSymbols))*D(i_lk,i_lk))).^2;        
            if mod(i_lk,1000)==0
                disp([int2str(i_lk/(obj.Nr.Subcarriers*obj.Nr.MCSymbols )*100) '%']);
            end 
        end      
        Frequency = (0:length(PowerSpectralDensity)-1)*1/(length(PowerSpectralDensity)*obj.PHY.dt);
        PowerSpectralDensity=PowerSpectralDensity/length(PowerSpectralDensity)^2/Frequency(2)^2;
        if nargout==0
            plot(Frequency,10*log10(PowerSpectralDensity));
            ylabel('Power Spectral Density (dB)');
            xlabel('Frequency (Hz)'); 
        end      
      end  
      
      function [PowerSpectralDensity,Frequency] = PlotPowerSpectralDensityUncorrelatedData(varargin)
        obj = varargin{1};        
        PowerSpectralDensity = zeros(obj.Nr.SamplesTotal,1);
        for i_lk = 1:obj.Nr.Subcarriers   
            V = zeros(obj.Nr.Subcarriers,obj.Nr.MCSymbols);
            V(i_lk,round(obj.Nr.MCSymbols/2))=1;
            PowerSpectralDensity = PowerSpectralDensity+abs(fft(obj.Modulation(V))).^2;        
        end   
        Frequency = (0:length(PowerSpectralDensity)-1)*1/(length(PowerSpectralDensity)*obj.PHY.dt);
        PowerSpectralDensity=obj.Nr.MCSymbols*PowerSpectralDensity/length(PowerSpectralDensity)^2/Frequency(2)^2;
        if nargout==0
            plot(Frequency,10*log10(PowerSpectralDensity));
            ylabel('Power Spectral Density (dB)');
            xlabel('Frequency (Hz)'); 
        end      
      end 

      %% SIR and Noise power  
      function Pn = GetSymbolNoisePower(varargin)
          obj = varargin{1};
          Pn_time = varargin{2};
          
          Pn = (Pn_time*(obj.Nr.Subcarriers*obj.PHY.SubcarrierSpacing/(obj.PHY.SamplingRate)));
      end 
      function [SignalPower,InterferencePower] = GetSignalAndInterferencePowerQAM(varargin)
          obj = varargin{1};
          VectorizedChannelCorrelationMatrix = varargin{2};
          DataSymbolCorrelationMatrix = varargin{3};           
          TimeSamplesOffset = varargin{4};
          if numel(varargin)==6          
            SubcarrierPosition = varargin{5};
            FBMCSymbolPosition = varargin{6};
          else
            SubcarrierPosition = 1:obj.Nr.Subcarriers;
            FBMCSymbolPosition = 1:obj.Nr.MCSymbols;
          end
          
          TXMatrix = obj.GetTXMatrix;
          RXMatrix = obj.GetRXMatrix;
          RXMatrix = [zeros(size(RXMatrix,1),TimeSamplesOffset),RXMatrix(:,1:end-TimeSamplesOffset)]; % Time offset compensation
          
          for i_l = 1:length(SubcarrierPosition)
          for i_k = 1:length(FBMCSymbolPosition)
             l = SubcarrierPosition(i_l);
             k = FBMCSymbolPosition(i_k);
             Index = l+(k-1)*obj.Nr.Subcarriers;
         
             %TempOld = kron(TXMatrix.',RXMatrix(Index,:))*VectorizedChannelCorrelationMatrix*kron(TXMatrix.',RXMatrix(Index,:))';
             % Much more efficient
             RXVectorRep = kron(sparse(eye(length(RXMatrix(Index,:)))),RXMatrix(Index,:)');
             Temp = RXVectorRep'*VectorizedChannelCorrelationMatrix*RXVectorRep;
             CorrMatrixNoData = TXMatrix.'*Temp*conj(TXMatrix);
             
             SignalDataSymbolCorrelationMatrix = zeros(size(DataSymbolCorrelationMatrix));
             SignalDataSymbolCorrelationMatrix(Index,Index) = DataSymbolCorrelationMatrix(Index,Index);
             InterferenceDataSymbolCorrelationMatrix = DataSymbolCorrelationMatrix;
             InterferenceDataSymbolCorrelationMatrix(Index,:)=0;
             InterferenceDataSymbolCorrelationMatrix(:,Index)=0;
            
             SignalPower(i_l,i_k) = abs(sum(sum(CorrMatrixNoData.*SignalDataSymbolCorrelationMatrix)));
             InterferencePower(i_l,i_k) = abs(sum(sum(CorrMatrixNoData.*InterferenceDataSymbolCorrelationMatrix)));
          end
          end
      end
   end 
end




%% Old stuff
% Super stupid
%          tic
%          RXVector = RXMatrix(Index,:);
%          LengthRXVector = length(RXVector);
%          Temp2 = nan(LengthRXVector,LengthRXVector);
%          for i_Row = 1:length(RXVector)
%          for i_Column = 1:length(RXVector)
%              IndexRow = (1:LengthRXVector)+(i_Row-1)*LengthRXVector;
%              IndexColumn = (1:LengthRXVector)+(i_Column-1)*LengthRXVector;
%              Temp2(i_Row,i_Column) = RXVector*VectorizedChannelCorrelationMatrix(IndexRow,IndexColumn)*RXVector';          
%          end
%          end
%          toc


% Not as stupid bust still bad
%              NumberEigenvalues = 6;
%              [U,S] = eigs(VectorizedChannelCorrelationMatrix,NumberEigenvalues);
%              while S(end,end)/sum(sum(S))>10^-3
%                NumberEigenvalues = NumberEigenvalues*2;
%               [U,S] = eigs(VectorizedChannelCorrelationMatrix,NumberEigenvalues);
%              end
%              SVDHelp = U*sqrt(S);
%              TempHelp = nan(size(TXMatrix,2),size(SVDHelp,2));
%              for i_SVDHelp = 1:size(SVDHelp,2)
%                  TempHelp(:,i_SVDHelp) = reshape((RXMatrix(Index,:)*reshape(SVDHelp(:,i_SVDHelp),size(TXMatrix,1),size(TXMatrix,1))*TXMatrix),[],1);
%              end
%              Temp = TempHelp*TempHelp';
