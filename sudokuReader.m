function sudokuReader(I,showSolution,templates,cornerExtractionMethoddigit,showProgress)

%% Pre-processing
I=rgb2gray(I);
SDK=double(I);

% Adjust contrast
SDK=imadjust(SDK.255).255;

% Threshold for binary image
SDK = SDK  mean2(SDK)  2;

% Remove small objects
SDK = ones(size(SDK))-SDK;
SDK_small_objects = SDK;
SDK = bwareaopen(SDK,floor(size(SDK,1)size(SDK,2)256));
SDK_small_objects = SDK_small_objects-SDK;

showImages(SDK,I,showProgress);

[h,w] = size(SDK);
if strcmp(cornerExtractionMethoddigit,'auto')
    % Extract corners (automatically)
    % Find a list of non-zero indices in SDK
    pList = find(SDK);
    % Extract the row-col indices of the above list
    [rowList colList] = ind2sub(size(SDK), pList);

    % Initialize all the corners to be the first point in the list
    tlRow = rowList(1); tlCol = colList(1);
    trRow = rowList(1); trCol = colList(1);
    brRow = rowList(1); brCol = colList(1);
    blRow = rowList(1); blCol = colList(1);

    % Go over all points and find the minimum vertical+horizontal
    % distance from closest image corners
    for i=1size(rowList)
        % Check for top left corner
        if rowList(i)+colList(i)tlRow+tlCol
            tlRow=rowList(i);
            tlCol=colList(i);
        end
        % Check for top right corner
        if rowList(i)+(w-colList(i))trRow+(w-trCol)
            trRow=rowList(i);
            trCol=colList(i);
        end
        % Check for bottom right corner
        if (h-rowList(i))+(w-colList(i))(h-brRow)+(w-brCol)
            brRow=rowList(i);
            brCol=colList(i);
        end
        % Check for bottom left corner
        if (h-rowList(i))+colList(i)(h-blRow)+blCol
            blRow=rowList(i);
            blCol=colList(i);
        end
    end

    % Plot the resulting corners
    subplot(1,2,1);
    plot(tlCol,tlRow,'b','linewidth',8);
    plot(trCol,trRow,'g','linewidth',8);
    plot(brCol,brRow,'r','linewidth',8);
    plot(blCol,blRow,'y','linewidth',8);

    % Plot the resulting corners
    subplot(1,2,2);
    plot(tlCol,tlRow,'b','linewidth',8);
    plot(trCol,trRow,'g','linewidth',8);
    plot(brCol,brRow,'r','linewidth',8);
    plot(blCol,blRow,'y','linewidth',8);

else
    % Extract corners (manually)
    title('Identify Top Left corner');
    [tlCol tlRow]=ginput(1);
    plot(tlCol,tlRow,'b','linewidth',8);

    title('Identify Top Right corner');
    [trCol trRow]=ginput(1);
    plot(trCol,trRow,'g','linewidth',8)

    title('Identify Bottom Right corner');
    [brCol brRow]=ginput(1);
    plot(brCol,brRow,'r','linewidth',8); 

    title('Identify Bottom Left corner');
    [blCol blRow]=ginput(1);
    plot(blCol,blRow,'y','linewidth',8)

end

if showProgress
    pause(0.5);
end

%% Calculate Homography and apply transformation to image

% Result list of 4 points
p = [tlCol, tlRow; trCol, trRow; brCol, brRow; blCol, blRow];

% second set of points is just a square
p2 = [0, 0; 0, w; h, w; h, 0];

A = zeros(8,8);
b=  zeros(8,1);
for i=14
    A(((i-1)2)+1, ) = [p(i,1) p(i,2) 1 0 0 0 -p(i,1)p2(i,1) -p(i,2)p2(i,1)];
    A(((i-1)2)+2, ) = [0 0 0 p(i,1) p(i,2) 1 -p(i,1)p2(i,2) -p(i,2)p2(i,2)];
    b(((i-1)2)+1) = p2(i,1);
    b(((i-1)2)+2) = p2(i,2);
end
h_vec = Ab;
H = [h_vec(13)'; h_vec(46)'; [h_vec(78)', 1]];

TFORM = maketform('projective', H');
[SDK XDATA YDATA] = imtransform(SDK,TFORM);
SDK = SDK';
SDK_small_objects = imtransform(SDK_small_objects,TFORM);
SDK_small_objects = SDK_small_objects';
I = imtransform(I,TFORM);
I = I';
ROWDATA=XDATA;
COLDATA=YDATA;

p(,1)=p2(,1)-ROWDATA(1);
p(,2)=p2(,2)-COLDATA(1);

tlRow = floor(max(1,p(1,1))); tlCol = floor(max(1,p(1,2)));
trRow = floor(max(1,p(2,1))); trCol = floor(min(size(SDK,2),p(2,2)));
brRow = floor(min(size(SDK,1),p(3,1))); brCol = floor(min(size(SDK,2),p(3,2)));
blRow = floor(min(size(SDK,1),p(4,1))); blCol = floor(max(1,p(4,2)));

showImages(SDK,I,false);

% Plot the new corners
subplot(1,2,1);
plot(tlCol,tlRow,'b','linewidth',8);
plot(trCol,trRow,'g','linewidth',8);
plot(brCol,brRow,'r','linewidth',8);
plot(blCol,blRow,'y','linewidth',8);

% Plot the new corners
subplot(1,2,2);
plot(tlCol,tlRow,'b','linewidth',8);
plot(trCol,trRow,'g','linewidth',8);
plot(brCol,brRow,'r','linewidth',8);
plot(blCol,blRow,'y','linewidth',8);

if showProgress
    pause(0.5);
end

%% Crop using new corners

SDK = SDK(tlRowbrRow,tlColbrCol);
SDK_small_objects = SDK_small_objects(tlRowbrRow,tlColbrCol);
I = I(tlRowbrRow,tlColbrCol);
showImages(SDK_small_objects,I,showProgress);

%% Divide into squares and identify numbers
[h,w] = size(SDK);
sqr_h = floor(h9);
sqr_w = floor(w9);
Board = zeros(9,9);

if sqr_h50
    filter_size = 2;
elseif sqr_h100
    filter_size = 3;
elseif sqr_h200
    filter_size = 4;
else
    filter_size = 5;
end

for row=19
    for col=19
        % Extract a square
        SQR = SDK_small_objects((row-1)sqr_h+1(row)sqr_h,(col-1)sqr_w+1(col)sqr_w);
        
        % Find and keep the biggest connected component
        CC = bwconncomp(SQR);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idx] = max(numPixels);
        SQR = zeros(size(SQR));
        SQR(CC.PixelIdxList{idx}) = 1;
        
        % Smoothen image and threshold
        SQR=conv2(double(SQR),fspecial('average',filter_size),'valid');
        SQR = SQR  0.9;
        
        % Remove anything that is not at least half the height of a square
        SQR = bwareaopen(SQR,floor(size(SQR,1)2));
        
        % Check if square contains a digit
        SQR_small = SQR(5size(SQR,1)-4,5size(SQR,2)-4);
        if size(find(SQR_small),1)==0
            % Empty square
            Board(row,col)=0;
        else
            % Extract digit, resize and compare
            [r,c] = find(SQR);
            n1 = SQR(min(r)max(r),min(c)max(c));  
            img_r = imresize(n1,[42 24]);
            Board(row,col) = read_letter(img_r,10,templates);
        end
    end
end

%% Display puzzle and solution
sqr_h = 64;
sqr_w = 64;
SDK = zeros(9sqr_h+8,9sqr_w+8);
SQR_template = zeros([sqr_h sqr_w]);

% Create square boundaries
SQR_template(1end,1) = 1;
SQR_template(1end,end) = 1;
SQR_template(1,1end) = 1;
SQR_template(end,1end) = 1;

% Solve the puzzle
if showSolution
    SudokuSolution = sudokuSolve(Board);
end

% Put numbers on board
for row=19
    for col=19
        SQR = SQR_template;
        if (showSolution) && Board(row,col)==0
            SQR(12sqr_h-11,21sqr_w-20) = templates{SudokuSolution(row,col)+1} . 3;
        else
            SQR(12sqr_h-11,21sqr_w-20) = templates{Board(row,col)+1};
        end
        
        % Create inside boundaries
        if mod(col,3)==0
            SQR(1end,end-2end-1) = 1;
        elseif mod(col,3)==1
            SQR(1end,23) = 1;
        end
        if mod(row,3)==0
            SQR(end-2end-1,1end) = 1;
        elseif mod(row,3)==1
            SQR(23,1end) = 1;
        end
        
        SDK(((row-1)sqr_h)+5rowsqr_h+4,((col-1)sqr_w)+5colsqr_w+4)=SQR;
    end
end

% Create board boundary
SDK(1end,14) = 1;
SDK(1end,end-3end) = 1;
SDK(14,1end) = 1;
SDK(end-3end,1end) = 1;

SDK = 1-SDK;
showImages(SDK,I,false);

end


% Simple Recursive Sudoku Solver 
function A=sudokuSolve(A) 
    [yy xx]=find(A==0); % find empty cells
    if isempty(xx)
        return 
    end
    x=xx(1); 
    y=yy(1); 
    for i=19 
        % try 1 to 9 
        y1=1+3floor((y-1)3); 
        % find 3x3 block 
        x1=1+3floor((x-1)3); 
        if ~( any(A(y, )==i)  any(A(,x)==i)  any(any(A(y1y1+2,x1x1+2)==i)) ) % if valid number 
            tmp=A; 
            tmp(y,x)=i; 
            tmp=sudokuSolve(tmp); % solve with this 
            if all(all(tmp)) 
                A=tmp;  % solution 
                return; % done 
            end
        end
    end
end

function digit=read_letter(imagn,num_letras,templates)
% Computes the correlation between template and input image.
% and its output is the best fit digit.
% Size of 'imagn' must be 42 x 24 pixels
    comp=[ ];
    for n=1num_letras
        sem=corr2(templates{1,n},imagn);
        comp=[comp sem];
    end
    digit=find(comp==max(comp))-1;
end

function showImages(SDK,I,showProgress)
    subplot(1,2,1); hold off;
    imshow(SDK); hold on;
    subplot(1,2,2); hold off;
    imshow(I); hold on;
    if showProgress
        pause(0.5);
    end
end