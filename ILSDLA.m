%% Iterative least squares (ILS) dictionary learning algorithm (DLA) with rapid OMP for representation
% Authored by Ayan Chatterjee (ayan@outlook.com) on 1 January 2020.

function [dictionary, sparseRep, t] = ILSDLA(Image, numAtoms, maxIter, maxAtomPerPixel)

%% pre-process
[h, w, b] = size(Image);
Image = reshape(single(Image), [h*w, b])'; % single datatype to save space
dictionary = single(randn(b, numAtoms));
dictionary = dictionary ./ sum(dictionary.^2).^0.5; % normalise the initial basis
t = zeros(1, maxIter);

% run the ILS-DLA
for iter = 1:maxIter
    tic;
    sparseRep = RapidOMP(Image, dictionary, maxAtomPerPixel);
    dictionary = Image*sparseRep'*pinv(sparseRep*sparseRep');
    dictionary = dictionary ./ sum(dictionary.^2).^0.5; % normalize the dictionary
    residue = Image - dictionary * sparseRep;
    t(iter) = toc; % captures time taken per training iteration
    disp(strcat("Iter: ", num2str(iter), ", mean residue: ", num2str(mean(sum(abs(residue)))), ", time taken: ", num2str(t(iter))));
end

sparseRep = reshape(sparseRep', [h, w, numAtoms]);

end