function [ output_args ] = PSC_LoadExcel(pathname, filename)
%% Function to read Excel file for PSC detection pipeline
    %Written by CRW, Sept 9 2018
        %last updated: Sept 10 2018


    evalin('base', 'global PSCTableRaw ')
    global PSCTableRaw 

 	clear PSCTableRaw
	evalin('base', 'global PSCTableRaw')
    global PSCTableRaw  
    
    %[filename,pathname,~] = uigetfile({'*.xlsx', '*xls*'}, 'Select excel file with annotations');
    
    if isequal(filename,0)
        disp('User selected Cancel')
        return
    else
        disp(['User selected ', fullfile(pathname, filename)])
        [PSCTableNum, PSCTableTxt, PSCTableRaw] = xlsread(fullfile(pathname, filename));
    end
    

    aa=cellfun(@isnan, PSCTableRaw(1,:), 'un',0);
    lastCol=find(cell2mat(cellfun(@all, aa, 'un',0)));
    if isempty(lastCol)
        lastCol=size(PSCTableRaw,2);
    else
        lastCol=lastCol(1)-1;
    end
    
    aa=cellfun(@isnan, PSCTableRaw(:,1), 'un',0);
    lastRow=find(cell2mat(cellfun(@all, aa, 'un',0)));
    if isempty(lastRow)
        lastRow=size(PSCTableRaw,1);
    else
        lastRow=lastRow(1)-1;
    end   
    
    disp(['Data loaded for ' num2str(lastRow-1) ' cells']);
    disp(' ');
    clear PSCTableNum PSCTableTxt
end 