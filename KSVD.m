%% K-SVD dictionary learning algorithm with rapid OMP for representation
% Authored by Ayan Chatterjee (ayan@outlook.com) on 24 December 2019. K-SVD algorithm was published with IEEE TSP journal (doi:10.1109/TSP.2006.881199).

function [dictionary, sparseRep, t] = KSVD(Image, numAtoms, maxIter, maxAtomPerPixel)

%% pre-process
[h, w, b] = size(Image);
Image = reshape(single(Image), [h*w, b])'; % single datatype to save space
dictionary = single(randn(b, numAtoms));
dictionary = dictionary ./ sum(dictionary.^2).^0.5; % normalise the initial basis
t = zeros(1, maxIter);

%% run K-SVD algorithm
for iter = 1:maxIter
    tic;
    sparseRep = RapidOMP(Image, dictionary, maxAtomPerPixel);
    residue = Image - dictionary * sparseRep;
    for atomPos = 1:numAtoms
        ImageSubsetPos = find(sparseRep(atomPos, :) > 0);
        if(numel(ImageSubsetPos) == 0)
            continue; % this atom is not contributing to any pixel
        end
        ImageSubset = residue(:, ImageSubsetPos) + dictionary(:, atomPos)*sparseRep(atomPos, ImageSubsetPos);
        [U, S, V] = svds(double(ImageSubset), 1, 'L'); % svds accepts double datatype
        dictionary(:, atomPos) = single(U); % U is already normalized
        sparseRep(atomPos, ImageSubsetPos) = S*V';
        residue(:, ImageSubsetPos) = ImageSubset - dictionary(:, atomPos)*sparseRep(atomPos, ImageSubsetPos);
    end
    t(iter) = toc; % captures time taken per training iteration
    disp(strcat("Iter: ", num2str(iter), ", mean residue: ", num2str(mean(sum(abs(residue)))), ", time taken: ", num2str(t(iter))));
end

sparseRep = reshape(sparseRep', [h, w, numAtoms]);

end