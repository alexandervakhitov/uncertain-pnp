function ShowInputImage(image1,v2d1,pointr,mask,pointp)
im1 = imread(image1);

% Show a figure with lines joining the accepted matches.
figure;
%imshow(im1);
imshow(im1,'border','tight','initialmagnification','fit');

hold on;

%show vertex projection 
if nargin > 1
    hold on; line([v2d1(1,1:4) v2d1(1,1)],[v2d1(2,1:4) v2d1(2,1)],'color','b','LineWidth',4);

    %for 3D box
    if size(v2d1,2) > 4
        hold on; line([v2d1(1,5:8) v2d1(1,5)],[v2d1(2,5:8) v2d1(2,5)],'color','b','LineWidth',4);
        hold on; line([v2d1(1,1) v2d1(1,5) v2d1(1,8) v2d1(1,4) v2d1(1,1)],[v2d1(2,1) v2d1(2,5) v2d1(2,8) v2d1(2,4) v2d1(2,1)],'color','b','LineWidth',4);
        hold on; line([v2d1(1,2) v2d1(1,6) v2d1(1,7) v2d1(1,3) v2d1(1,2)],[v2d1(2,2) v2d1(2,6) v2d1(2,7) v2d1(2,3) v2d1(2,2)],'color','b','LineWidth',4);
    end
    
    hold on; plot(v2d1(1,1),v2d1(2,1),'o','MarkerSize',15,'MarkerFaceColor','m');
end

%show the inliers and outliers
if nargin > 2 
    for i = 1:size(pointr,2)
        if mask(i) > 0
            hold on; plot(pointr(1,i),pointr(2,i),'c+','MarkerSize',10); 
        else
            hold on; plot(pointr(1,i),pointr(2,i),'r+','MarkerSize',10);
        end
    end
end

%show the projection of inliers
if nargin > 4 
    for i = 1:size(pointp,2)
        hold on; plot(pointp(1,i),pointp(2,i),'go','MarkerSize',10);
    end
end
axis off;

hold off;
end

