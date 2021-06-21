function ShowInlierOutlier(point1,point2,image1,image2,y,v2d1,v2d2)
im1 = imread(image1);
im2 = imread(image2);

im3 = appendimages(im1,im2);

% Show a figure with lines joining the accepted matches.
figure('Position', [100 100 size(im3,2) size(im3,1)]);
colormap('gray');
imagesc(im3);
cols1 = size(im1,2);
hold on;

%show vertex projection 
if nargin > 5
    hold on; plot(v2d1(1,1),v2d1(2,1),'o','MarkerSize',15,'MarkerFaceColor','c');

    hold on; line([v2d1(1,1:4) v2d1(1,1)],[v2d1(2,1:4) v2d1(2,1)],'color','b','LineWidth',4);

    %for 3D box
    if size(v2d2,2) > 4
        hold on; line([v2d1(1,5:8) v2d1(1,5)],[v2d1(2,5:8) v2d1(2,5)],'color','b','LineWidth',4);
        hold on; line([v2d1(1,1) v2d1(1,5) v2d1(1,8) v2d1(1,4) v2d1(1,1)],[v2d1(2,1) v2d1(2,5) v2d1(2,8) v2d1(2,4) v2d1(2,1)],'color','b','LineWidth',4);
        hold on; line([v2d1(1,2) v2d1(1,6) v2d1(1,7) v2d1(1,3) v2d1(1,2)],[v2d1(2,2) v2d1(2,6) v2d1(2,7) v2d1(2,3) v2d1(2,2)],'color','b','LineWidth',4);
    end
    
    v2d2(1,:) = v2d2(1,:) + cols1; %in new coordinate system
    
    hold on; plot(v2d2(1,1),v2d2(2,1),'o','MarkerSize',15,'MarkerFaceColor','c');

    hold on; line([v2d2(1,1:4) v2d2(1,1)],[v2d2(2,1:4) v2d2(2,1)],'color','b','LineWidth',4);
      
    if size(v2d2,2) > 4
        hold on; line([v2d2(1,5:8) v2d2(1,5)],[v2d2(2,5:8) v2d2(2,5)],'color','b','LineWidth',4);
        hold on; line([v2d2(1,1) v2d2(1,5) v2d2(1,8) v2d2(1,4) v2d2(1,1)],[v2d2(2,1) v2d2(2,5) v2d2(2,8) v2d2(2,4) v2d2(2,1)],'color','b','LineWidth',4);
        hold on; line([v2d2(1,2) v2d2(1,6) v2d2(1,7) v2d2(1,3) v2d2(1,2)],[v2d2(2,2) v2d2(2,6) v2d2(2,7) v2d2(2,3) v2d2(2,2)],'color','b','LineWidth',4);
    end
end


for i = 1: size(point1,2)
% figure('Position', [100 100 size(im3,2) size(im3,1)]);
% colormap('gray');
% imagesc(im3);
% hold on;
% cols1 = size(im1,2);
if y(i) == 1
        line([point1(1,i) point2(1,i)+cols1], ...
              [point1(2,i) point2(2,i)], 'Color', 'g','LineWidth',2);
else

     line([point1(1,i) point2(1,i)+cols1], ...
             [point1(2,i) point2(2,i)], 'Color', 'r','LineWidth',4, 'LineStyle','-.');

end
%        close all;

end
axis off;
hold off;

