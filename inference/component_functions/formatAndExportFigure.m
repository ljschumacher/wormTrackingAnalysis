function [] = formatAndExportFigure(handle,figurename,exportOptions)
set(handle,'PaperUnits','centimeters')
exportfig(handle,[figurename '.eps'],exportOptions)
system(['epstopdf ' figurename '.eps']);
system(['rm ' figurename '.eps']);
end