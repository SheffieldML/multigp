function m = lfmMeanCompute(meanFunction, X, varargin)

% LFMMEANCOMPUTE Give the output of the lfm mean function model for given X.
% FORMAT
% DESC gives the output of the lfm mean function model for a given input X.
% ARG model : structure specifying the model.
% ARG X : input location(s) for which output is to be computed.
% RETURN Y : output location(s) corresponding to given input
% locations.
%
% SEEALSO : lfmMeanCreate
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008


% MULTIGP

nlf = varargin{1};
startVal=1;
endVal=0;
for i =1:nlf,
    endVal = endVal + size(X{i}, 1);
    m(startVal:endVal, 1) =  zeros(length(X{i}),1);
    startVal = endVal+1;
end

for i = nlf+1:length(X),
    endVal = endVal + size(X{i}, 1);
    m(startVal:endVal, 1) = meanFunction.basal(i-nlf)/meanFunction.spring(i-nlf) * ...
        ones(length(X{i}),1);
    startVal = endVal+1;
end