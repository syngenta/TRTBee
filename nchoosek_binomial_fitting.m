function out = nchoosek_binomial_fitting(nVec, kVec)

for i = 1:length(nVec)
    out(i) = nchoosek(nVec(i), kVec(i));
end