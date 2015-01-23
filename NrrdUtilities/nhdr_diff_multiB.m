function [S, u, b, voxel, Sm, s0] = nhdr_diff_multiB(fn)

 %output
 % S - the DWI data
 % u - the gradient directions
 % b - the b-value
 % voxel - the voxel size
 % Sm - the space directions matrix
 
%  paths;
  
  dwi = loadNrrdStructureMultiB(fn);
  S = dwi.data;
  clear dwi.data;
  
  %u = dwi.gradients;
  M = dwi.measurementframe;
  bmax = dwi.bvalue;
  order = [1 2 3 4];
  %make adjustment of the spacedirections and measurement frame
  
  sd = dwi.spacedirections;
  
% %   %if space dimentions are flipped in either the x or y axis
% %   %don't know why z would still work !
% %   dsd = (diag(sd) < 0);
% %   if dsd(1) 
% %       S = flipdim(S,1);
% %   end
% %   if dsd(2)
% %       S = flipdim(S,2);
% %   end
  
  voxel = [norm(sd(:,1));norm(sd(:,2));norm(sd(:,3))]';
  R = sd./repmat(voxel,[3 1]);
  %b = sqrt(sum(dwi.gradients.^2,2)); %this is for phantom data, and MGH
  b = (sum(dwi.gradients.^2,2)); %for in-vivo data
  
  dwi.gradients(b>=0.05,:) = dwi.gradients(b>=0.05,:)./sqrt(b(b>=0.05,ones(1,3)));

  u = inv(R)*M*dwi.gradients';u=u';
  %b  = b.^2;
  
  %Sm is the space directions needed to transform the fiber tracts in the
  %slicer space
  if strcmp(dwi.space,'LPS') || strcmpi(dwi.space,'left-posterior-superior')
    RAS = [-1 0 0;0 -1 0;0 0 1];
    Sm = RAS * dwi.spacedirections;
  else
    Sm = dwi.spacedirections;
  end
  %Sm = inv(R)*M;
  
  ndirs = size(u,1);
  
  %find the storage format
  switch size(u,1)
      case size(S,1)
         order = [2 3 4 1];
      case size(S,2)
         order = [1 4 3 2];
      case size(S,3)
          order = [1 2 4 3];
  end
  S = permute(S,order);
  sz = size(S);
  Sm(4,:) = -0.5 * Sm * sz(1:3)'; %the space origin in Slicer
  
  
  % first find the number of baseline data
  idb=[];
  id=[];
  for j =1:ndirs
      %the baseline images
      if norm(u(j,:)) <= 0.05
          idb = [idb;j];
      else
          id = [id;j];
      end
  end
  
  if isempty(idb)
      fprintf('Data does not have a B0 image. \n');
       if max(b) < 1.1
        b = b.*bmax(id)';
       end
      return;
  end
  ndir = length(id);
  
 
  % separate out baseline
  s0 = S(:,:,:,idb);
  S  = S(:,:,:,id); % drop null-gradient slices
  
  % actual gradient directions used
  u = u(id,:);
  u = [u; -u]; % antipodal
  b = flat(b(id));
  if max(b) < 1.1
      b = b.*bmax(id)';
  end
  b=[b;b];
  % divide off baseline
  s0 = single(mean(s0, 4)); if nargout == 4, baseline = s0; end
  
  s0(s0==0) = 1; 
  %S = single(S) ./ s0(:,:,:,ones(1,ndir));

  S = single(S);
  for j=1:size(S,4)
      S(:,:,:,j) = S(:,:,:,j)./s0;
  end
  
end
