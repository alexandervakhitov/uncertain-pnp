function ShowReferenceImage(image1,v2d1)
im1 = imread(image1);

% Show a figure with lines joining the accepted matches.
figure;
%imshow(im1);
imshow(im1,'border','tight','initialmagnification','fit');

hold on;

%show vertex projection 
if nargin > 1
   
    hold on; line([v2d1(1,1:4) v2d1(1,1)],[v2d1(2,1:4) v2d1(2,1)],'color','c','LineWidth',4);

    %for 3D box
    if size(v2d1,2) > 4
        hold on; line([v2d1(1,5:8) v2d1(1,5)],[v2d1(2,5:8) v2d1(2,5)],'color','c','LineWidth',4);
        hold on; line([v2d1(1,1) v2d1(1,5) v2d1(1,8) v2d1(1,4) v2d1(1,1)],[v2d1(2,1) v2d1(2,5) v2d1(2,8) v2d1(2,4) v2d1(2,1)],'color','c','LineWidth',4);
        hold on; line([v2d1(1,2) v2d1(1,6) v2d1(1,7) v2d1(1,3) v2d1(1,2)],[v2d1(2,2) v2d1(2,6) v2d1(2,7) v2d1(2,3) v2d1(2,2)],'color','c','LineWidth',4);
    end
    
     hold on; plot(v2d1(1,1),v2d1(2,1),'o','MarkerSize',15,'MarkerFaceColor','m');
end
axis off;

%set(gca, 'Position', get(gca, 'OuterPosition') - ...
%     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
 
hold off;
end

