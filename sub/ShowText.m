function T = ShowText(str,x,y,fs,col)
% ShowText(str) shows text to gca on default position.
% ShowText(str,x) specifies x positions
% ShowText(str,x,y) specifies x and y positions
% ShowText(str,x,y,fs) specifies fontsize
% ShowText(str,x,y,fs,col) specifies color, color is triplet or character
% You can use default value like ShowText('test string',[],[],[],'r')

% default setting if no input argument
xlim1 = get(gca,'Xlim');
ylim1 = get(gca,'Ylim');
xPos  = xlim1(1) + diff(xlim1)*7/10;
yPos  = ylim1(1) + diff(ylim1)*3/5;
global defaultFontSize
if ~isempty(defaultFontSize)
    fontsize = defaultFontSize;
else
    fontsize = 20;
end
color = [0 0 0];

switch nargin
    case 2       
        if ~isempty(x)
            xPos = x;
        end
    case 3
        if ~isempty(x)
            xPos = x;
        end
        if ~isempty(y)
            yPos = y;
        end
    case 4
        if ~isempty(x)
            xPos = x;
        end
        if ~isempty(y)
            yPos = y;
        end
        if ~isempty(fs)
            fontsize = fs;
        end
    case 5
        if ~isempty(x)
            xPos = x;
        end
        if ~isempty(y)
            yPos = y;
        end
        if ~isempty(fs)
            fontsize = fs;
        end
        if ~isempty(col)
            color = col;
        end
end
    
    T = text(xPos, yPos, str, 'color', color);
    set(T, 'FontSize', fontsize);    
end
