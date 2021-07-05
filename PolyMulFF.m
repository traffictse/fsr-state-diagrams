function varargout=PolyMulFF(modulus,varargin)
if nargin==1
    error('Not enough inputs! At least 2!');
end
if nargin==2
        varargout{1}=varargin{1};
end
if nargin>2
    if (strcmp(varargin{nargin-2},'power'))
        pow=varargin{nargin-1};
        if(length(pow)<nargin-3)
            pow(length(pow)+1:nargin-3)=1;
        end
    else if (strcmp(varargin{nargin-1},'power'))
            error('Not enough inputs to be powers!');
        else
            pow(1:nargin-1)=1;
        end
    end
    product=1;
    for i=1:length(pow)
        for j=1:pow(i)
            product=gfconv(product,varargin{i},modulus);
        end
    end
    varargout{1}=product;
end
end