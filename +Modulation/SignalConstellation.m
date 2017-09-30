classdef SignalConstellation < handle 
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2016 by Institute of Telecommunications, TU Wien
% www.tc.tuwien.ac.at    

  properties (SetAccess = private)
       Method
       ModulationOrder
       BitMapping
       SymbolMapping
       Implementation
  end
  
  methods
  	function obj = SignalConstellation(varargin)
                
        obj.ModulationOrder = varargin{1};
        obj.Method = varargin{2};
        
        %% Gray coded bitmapping
        if strcmp(obj.Method,'QAM' )
            % Old method which uses model object
%             modulation_object=modem.qammod('M',obj.ModulationOrder,'SymbolOrder','gray');
            BitMappingAtom = [ones(sqrt(obj.ModulationOrder)/2,1);zeros(sqrt(obj.ModulationOrder)/2,1)];
            for i_temp = 2:log2(sqrt(obj.ModulationOrder))
                BinaryTemp = BitMappingAtom(1:2:end,i_temp-1);
                BitMappingAtom(:,i_temp) = [BinaryTemp;BinaryTemp(end:-1:1)];
            end
            IQ = 2*(1:sqrt(obj.ModulationOrder))-sqrt(obj.ModulationOrder)-1;
            [I_rep,Q_rep]=meshgrid(IQ,IQ);
            obj.SymbolMapping = I_rep(:)+1i*Q_rep(:);
            obj.SymbolMapping = obj.SymbolMapping/sqrt(mean(abs(obj.SymbolMapping).^2));
            obj.BitMapping = false(obj.ModulationOrder,log2(obj.ModulationOrder));
            for x_IQ = IQ
                obj.BitMapping(I_rep(:)==x_IQ,2:2:end) = BitMappingAtom;
                obj.BitMapping(Q_rep(:)==x_IQ,1:2:end) = BitMappingAtom;
            end
        elseif strcmp(obj.Method,'PAM')
%           modulation_object=modem.pammod('M',obj.ModulationOrder,'SymbolOrder','gray');        
            obj.BitMapping = [ones(obj.ModulationOrder/2,1);zeros(obj.ModulationOrder/2,1)];
            for i_temp = 2:log2(obj.ModulationOrder)
                BinaryTemp = obj.BitMapping(1:2:end,i_temp-1);
                obj.BitMapping(:,i_temp) = [BinaryTemp;BinaryTemp(end:-1:1)];
            end
            obj.SymbolMapping = (2*(1:obj.ModulationOrder)-obj.ModulationOrder-1).';
            obj.SymbolMapping = obj.SymbolMapping/sqrt(mean(abs(obj.SymbolMapping).^2));
        else
           error('Signal constellation method must be QAM or PAM!');
        end
        
        [~,SortOrder] = sort(bi2de(obj.BitMapping),'ascend');
        obj.SymbolMapping = obj.SymbolMapping(SortOrder);
        obj.BitMapping = obj.BitMapping(SortOrder,:);
        
%       	obj.SymbolMapping=modulation_object.modulate(0:obj.ModulationOrder-1);
%         Normalization=sqrt(mean(abs(obj.SymbolMapping).^2));
%         obj.BitMapping=de2bi(0:obj.ModulationOrder-1,log2(obj.ModulationOrder));
%         obj.SymbolMapping=obj.SymbolMapping(:)./Normalization;

%        For QAM modulation order >= 256, delaunayTriangulation performs
%        better than the method here 
%         obj.DT=delaunayTriangulation([real(obj.SymbolMapping(:)) imag(obj.SymbolMapping(:))]);

        obj.Implementation.DataSymbolsBitvalueOne = repmat(obj.SymbolMapping,1,log2(obj.ModulationOrder));
        obj.Implementation.DataSymbolsBitvalueOne=reshape(obj.Implementation.DataSymbolsBitvalueOne(logical(obj.BitMapping)),obj.ModulationOrder/2,[]);
        obj.Implementation.DataSymbolsBitvalueZero = repmat(obj.SymbolMapping,1,log2(obj.ModulationOrder));
        obj.Implementation.DataSymbolsBitvalueZero=reshape(obj.Implementation.DataSymbolsBitvalueZero(not(logical(obj.BitMapping))),obj.ModulationOrder/2,[]);       
        
    end   
    
    function DataSymbols = Bit2Symbol(varargin)
        obj = varargin{1};
        BinaryStream = varargin{2};
        DataSymbols = obj.SymbolMapping(bi2de(reshape(BinaryStream,log2(obj.ModulationOrder),[])')+1);            
    end

    function EstimatedBitStream = Symbol2Bit(varargin)
        obj = varargin{1};
        EstimatedDataSymbols = varargin{2};
        EstimatedDataSymbols = EstimatedDataSymbols(:);
        
        % only efficient for small modulation orders, i.e., <= 64
        % maybe rewrite it
        [~,b] =min(abs((repmat(EstimatedDataSymbols,1,obj.ModulationOrder)-repmat((obj.SymbolMapping).',size(EstimatedDataSymbols,1),1)).'));
        EstimatedBitStream = obj.BitMapping(b(:),:).';
        EstimatedBitStream = EstimatedBitStream(:);         
    end
    
    function QuantizedDataSymbols = SymbolQuantization(varargin)
        obj = varargin{1};
        EstimatedDataSymbols = varargin{2};
        EstimatedDataSymbols = EstimatedDataSymbols(:);
        
        % only efficient for small modulation orders, i.e., <= 64
        % maybe rewrite it
        [~,b] =min(abs((repmat(EstimatedDataSymbols,1,obj.ModulationOrder)-repmat((obj.SymbolMapping).',size(EstimatedDataSymbols,1),1)).'));
        QuantizedDataSymbols = obj.SymbolMapping(b(:),:).';
        QuantizedDataSymbols = QuantizedDataSymbols(:);
    end    
    
    function LLR = LLR_AWGN(varargin)
        % This method ...
        obj = varargin{1};
        ReceivedDataSymbols = varargin{2};
        Pn = varargin{3};
    
        if numel(Pn)>1
            PnRepeated = reshape(repmat(Pn.',log2(obj.ModulationOrder)*obj.ModulationOrder/2,1),obj.ModulationOrder/2,[]);
        else 
            PnRepeated = Pn;
        end
  
        ReceivedDataSymbolsRepeated = reshape(repmat(ReceivedDataSymbols.',log2(obj.ModulationOrder)*obj.ModulationOrder/2,1),obj.ModulationOrder/2,[]);
        DataSymbolsBitvalueOneRepeated = repmat(obj.Implementation.DataSymbolsBitvalueOne,1,length(ReceivedDataSymbols));
        DataSymbolsBitvalueZeroRepeated = repmat(obj.Implementation.DataSymbolsBitvalueZero,1,length(ReceivedDataSymbols));
        LLR = log(sum(exp(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueOneRepeated).^2./PnRepeated),1)./...
            sum(exp(-abs(ReceivedDataSymbolsRepeated-DataSymbolsBitvalueZeroRepeated).^2./PnRepeated),1)).';
        LLR(LLR==Inf)=10^10;
        LLR(LLR==-Inf)=-10^10;     
    end


     
  end
end
      