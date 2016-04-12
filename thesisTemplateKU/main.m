

% ear recognition demo
str{1}='back';str{2}='down';str{3}='front';str{4}='left';str{5}='right';str{6}='up';str{7}='zoom';
% num1=48;% num can not be equal to following numbers due to lack of data:
% num2=46;% 6,15,16,17,49,50,60 the maximal value is 106
% flag1=3;
% flag2=6;
% filename1=strcat('.\subset-1\',sprintf('%.3d',num1),'_',str{flag1},'_ear.jpg');
% filename2=strcat('.\subset-1\',sprintf('%.3d',num2),'_',str{flag2},'_ear.jpg');
num1=2;% number should be less than 125
num2=2;
flag1=3;
flag2=5;% flag number should be less than three
filename1=strcat('.\raw\',sprintf('%.3d',num1),'_',num2str(flag1),'.bmp');
filename2=strcat('.\raw\',sprintf('%.3d',num2),'_',num2str(flag2),'.bmp');
[im1, des1, loc1] = sift(filename1);
[im2, des2, loc2] = sift(filename2);
[num,matches] = match(im1, des1, loc1, im2, des2, loc2) ;
if length(matches)<10
    result=0;
    fprintf('Matching Failed!\n');
    title('Matching Failed!')
    return
end
x1=[loc1(matches(1,:),1:2),ones(size(matches,2),1)]';
x2=[loc2(matches(2,:),1:2),ones(size(matches,2),1)]';
[F,inliers]=RobustFundaMatrix(x1,x2);
threshold=0.8;
if length(inliers)/size(matches,2)>threshold
    result=1;
    fprintf('Matching Succeed!\n');
    title(['Matching Succeeded! ',' | Number of Matches are: ',num2str(length(matches)),' , Percentage is : ',num2str((length(inliers)/size(matches,2))*100),'%'])
else
    result=0;
    fprintf('Matching Failed!\n');
    title(['Matching Failed ! ','| Number of Mathces are:',num2str(length(matches)),' , Percentage is : ',num2str((length(inliers)/size(matches,2))*100),'%'])
end
