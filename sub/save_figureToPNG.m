function [] = save_figureToPNG( handle, filePathName, paper_DPI )

	% Determine the size of PDF output in terms of paper size (in inches)
	screen_DPI = get( 0, 'ScreenPixelsPerInch' );
	pos = getpixelposition( handle );
	imageSize = pos(3:4);						% in pixels
	paperSize = imageSize / screen_DPI * 2.54;	% in cm

	% Set the proper properties
	set( handle, 'PaperUnits', 'centimeters' );
	set( handle, 'PaperSize', paperSize );
%	set( handle, 'PaperPosition', [0 0 paperSize] );
 	set( handle, 'PaperPositionMode', 'auto' );
	
	% Print out to a PNG file (NOTE: Without '-r' option, default DPI for PNG format is 150.)
	if ( exist('paper_DPI','var') )
		print( handle, filePathName, '-dpng', ['-r' num2str(paper_DPI)] );
	else
		print( handle, filePathName, '-dpng', ['-r' num2str(screen_DPI)] );
	end
end