[image1,descrips,locs]=sift('ear1.pgm');
% You cannot change the name of the locs, that is why,it is being stored in
% loc1 for later use
locs1=locs;%loc1 is locs for first first image,same with loc2 and descrips1 and descrips 2.
descrips1=descrips;
[image2,descrips,locs]=sift('ear1.pgm');
locs2=locs;
descrips2=descrips;
[num,matches] = match(image1,descrips1,locs1,image2,descrips2,locs2);
if num<5
    fprintf('Match Failed')
    title('Match Failed')
    return
end
x1=[loc1(matches(1,:),1:2),ones(size(matches,2),1)]';
x2=[loc2(matches(2,:),1:2),ones(size(matches,2),1)]';
[F,inliers]=RobustFundaMatrix(x1,x2);
threshold=0.7;
if length(inliers)/size(matches,2)>threshold
    result=1;
    fprintf('Match');
    title('Match')
else
    result=0;
    fprintf('Failed');
    title('Failed');
end






