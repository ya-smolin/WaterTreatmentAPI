classdef View < handle
    
    properties
        hGUI;
        hModel;
    end
    
    properties(SetObservable = true)
        tableRows;%?need
        hDataAxes
    end
    
    methods
        
        function this = View(Model)
            %VIEW  a GUI representation of the signal model
            handles = View.initGUI();
            
            this.tableRows = this.formIsotermTable(Model, handles.tabOut);
            isotermsListener = event.proplistener(Model, findprop(Model,'isoterms'), 'PostSet',...
                @(o, e)this.update(o, e));
            setappdata(handles.fig, 'proplistener',isotermsListener);
            
            modelDataListener = event.proplistener(Model, findprop(Model,'data'), 'PostSet',...
                @(o, e)this.updateData());
            setappdata(handles.fig, 'proplistener',modelDataListener);
            
            this.hGUI = handles;
            this.hModel = Model;
        end
        
        function update(this, o, e)
            isotermId = this.hModel.lastIsoInd;
            if(isotermId == -1)
                this.formIsotermTable(this.hModel, this.hGUI.tabOut);
                return;
            end
            isoterm = this.hModel.isoterms{isotermId};
            if(isempty(isoterm) || isempty(isoterm.isotermResult))
                isoterm = this.hModel.isotermTypes{isotermId};
            end
            
            %clear plot
            if(~isempty(this.tableRows{isotermId}))
                this.tableRows{isotermId}.clearPlots();
            end
            
            confInt = str2double(get(this.hGUI.edConfInt, 'String')) / 100;
            tableRow = IsotermTableRow(isoterm, confInt);
            %update table
            table = this.hGUI.tabOut;
            tableRowsData = get(table, 'Data');
            tableRowsData(isotermId, :) = tableRow.data;
            set(table, 'Data', tableRowsData);
            
            %update plot
            if(~isempty(isoterm) && isprop(isoterm, 'isotermResult'))
                isChecked =  get(this.hGUI.cbConfInt, 'Value');
                if(isChecked)
                    tableRow.hConfIntAxes = this.plotConfInt(isoterm);
                end
                tableRow.axisHandles = plot(isoterm.isotermResult);
                tableRow.initHoverLegend();
                legend('off');
            end
            this.tableRows{isotermId} = tableRow;
        end
        
        function hConfPlot = plotConfInt(this, isoterm)
            xConfInt = linspace(0, max(this.hModel.data(:,1)), 400);
            confidenceLevel = str2double(get(this.hGUI.edConfInt,'String')) / 100;
            yConfInt = predint(isoterm.isotermResult, xConfInt, confidenceLevel);
            yConfInt(yConfInt<0) = 0;
            hConfPlot = plot(xConfInt, yConfInt, '--c');
        end
        
        function updateData(this)
            data = this.hModel.data;
            delete(this.hDataAxes);
            cla;
            handles = this.hGUI;
            axes(handles.axes);
            data(end, :) = [0, 0];
            Cr = data(:, 1);
            Ar = data(:, 2);
            this.hDataAxes = plot(Cr, Ar, 'greeno', 'LineWidth', 3);
            for i = 1:this.hModel.size
                this.tableRows{i}.axisHandles = -1;
                this.tableRows{i}.hConfIntAxes = -ones(1, 2);
            end
            
            hold on;
        end
        
        function isotermIdList = getCheckedRows(this)
            tableOut = this.hGUI.tabOut;
            dataOut = get(tableOut, 'Data');
            isShownColumn = dataOut(:, IsotermTableRow.columnShow);
            isotermIdList = find(cell2mat(isShownColumn) == 1)';
        end
        
        function isotermIdList = getAllCalculatedIsoterms(this)
            isoterms = this.hModel.isoterms;
            isotermIdList = [];
            for i = 1:Model.size
                isot = isoterms{i};
                if(~isempty(isot) && ~isempty(isot.isotermResult))
                    isotermIdList(end+1) = i;
                end;
            end
        end
    end
    
    methods(Static)
        function tableRows = formIsotermTable(M, isotermTable)
            N = length(M.isoterms);
            tableRows = cell(N, 1);
            tableRowsData = cell(N, IsotermTableRow.size);
            for i = 1:N
                if isempty(M.isoterms{i})
                    tableRows{i} = IsotermTableRow(M.isotermTypes{i});
                else
                    tableRows{i} = IsotermTableRow(M.isoterms{i});
                end
                tableRowsData(i, :) = tableRows{i}.data;
            end
            set(isotermTable, 'Data', tableRowsData,...
                'ColumnName', IsotermTableRow.columnName, 'ColumnFormat', IsotermTableRow.columnFormat,...
                'ColumnEditable', IsotermTableRow.columnEditable);
            
            jscrollpane = findjobj(isotermTable);
            jtable = jscrollpane.getViewport.getView;
            
            % Now turn the JIDE sorting on
            jtable.setSortable(true);		% or: set(jtable,'Sortable','on');
            jtable.setAutoResort(true);
            jtable.setPreserveSelectionsAfterSorting(true);
            jtable.setMultiColumnSortable(true);
            
            
            View.fitTableWidth(isotermTable);
           
        end
        
        function fitTableWidth(table)
%             jscrollpane = findjobj(table);
%             jtable = jscrollpane.getViewport.getView;
%             jtable.setAutoResizeMode(jtable.AUTO_RESIZE_SUBSEQUENT_COLUMNS);
            tableRowsData = get(table, 'Data');
            columnWidth = max(cellfun('length', tableRowsData));
            columnWidth(columnWidth < 6) = 0;
            columnWidth = columnWidth * 7;
            columnWidthCell = num2cell(columnWidth);
            for i = 1:length(columnWidthCell)
                if(columnWidthCell{i} == 0)
                    columnWidthCell{i} = 'auto';
                end
            end
            set(table,'ColumnWidth',columnWidthCell);
        end
        
        function handles = initGUI()
            % load FIG file (its really a MAT-file)
            hFig = hgload('View.fig');
            handles = struct('fig',hFig);
            tags = {'axes', 'cbConfInt', 'edConfInt', 'butRecalc', 'tabIn',...
                'tabOut', 'butPlusRow', 'butMinusRow', 'butLoadData',...
                'menuTabIn', 'contextPaste', 'edP1', 'tvP1', 'cbP1', 'barSave', 'barDataCursor', 'barZoomIn', 'barZoomOut'};
            % extract handles to GUI components
            for tag = tags
                tag = char(tag);
                hTag = findobj(hFig, 'tag', tag);
                handles.(tag) = hTag;
            end
            
            set(handles.tabIn, 'uicontextmenu', handles.menuTabIn);
        end
        
    end
end
