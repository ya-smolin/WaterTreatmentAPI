classdef Controller < handle
    properties
        model;
        view;
        data = [];
    end
    methods
        function this = Controller()
            model = Model();
            view = View(model);
            set(view.hGUI.contextPaste, 'Callback',{@Controller.copyExcel,view.hGUI.tabIn});
            set(view.hGUI.butPlusRow, 'Callback',{@Controller.onClickPlusRow,view.hGUI.tabIn});
            set(view.hGUI.butMinusRow, 'Callback',{@Controller.onClickMinusRow,view.hGUI.tabIn});
            set(view.hGUI.butLoadData, 'Callback',@this.onClickLoadData);
            set(view.hGUI.tabOut, 'CellEditCallback', @this.onTableEdit);
            set(view.hGUI.butRecalc, 'Callback', @this.onClickRecalculate);
            set(view.hGUI.cbP1, 'Callback', @this.onCheckP1);
            set(view.hGUI.cbConfInt, 'Callback', @this.onCheckConfInt);
            set(view.hGUI.edConfInt, 'Callback', @this.onFinishConfInt);
            %2,4DNP
            Cr = [
                0.5
                1
                2
                4
                8
                18
                46
                50
                80
                90
                ];
            
            Ar = [
                83.25
                99.8
                124.5
                165.3333333
                246
                337.4
                385.9
                450
                504
                546.6666667
                ];
            data = [Cr Ar];
            set(view.hGUI.tabIn, 'data', data);
            this.model = model;
            this.view = view;
        end
        function onFinishConfInt(this, edConfInt, event)
            confAxes = this.view.hConfIntAxes;
            isotermsId = this.view.getCheckedRows();
            for id = isotermsId
                hConfInt = confAxes(id, :);
                if(ishghandle(hConfInt))
                    delete(hConfInt);
                end
                this.view.hConfIntAxes(id, :) = this.view.plotConfInt(this.model.isoterms{id});
            end
        end
        
        function onCheckConfInt(this, cbConfInt, event)
            isChecked =  get(this.view.hGUI.cbConfInt, 'Value');
            chekedIsotermsID = this.view.getCheckedRows();
            for isotermId = chekedIsotermsID
                if(isChecked)
                    isoterm = this.model.isoterms{isotermId};
                    if ~isempty(isoterm)
                        hConfPlot = this.view.plotConfInt(isoterm);
                        this.view.hConfIntAxes(isotermId, :) = hConfPlot;
                    end
                else
                    hConfLine = this.view.hConfIntAxes(isotermId, :);
                    if(ishghandle(hConfLine)) 
                        delete(hConfLine);
                    end;
                end
            end
        end
        
        function onCheckP1(this, cbP1, event)
            isChecked =  get(cbP1, 'Value');
            hEdP1 = this.view.hGUI.edP1;
            if(isChecked)
               set(hEdP1, 'Visible', 'off');
            else
               set(hEdP1, 'Visible', 'on');
            end
        end
        
        function onClickRecalculate(this, butRecalc, event)
            isChecked =  get(this.view.hGUI.cbP1, 'Value');
            if(isChecked)
                this.model.parameters = [];
            else
                this.model.parameters = str2num(get(this.view.hGUI.edP1, 'String'));
            end
            isotermsIdList = this.view.getCheckedRows();
            tableIn = this.view.hGUI.tabIn;
            dataIn = get(tableIn, 'Data');
            if(~isequal(this.data, dataIn))
                this.model.data = dataIn;
                this.model.clear();
            end
            this.data = dataIn;
            
            table = this.view.hGUI.tabOut;
            tableRowsData = get(table, 'Data');
            for i = isotermsIdList
                tableRowsData{i, IsotermTableRow.columnShow} = true;
            end
            set(table, 'Data', tableRowsData);
            
            this.model.calculate(isotermsIdList);
        end

        function onTableEdit(this, table, eventData)
            row = eventData.Indices(1);
            col = eventData.Indices(2);
            tableOut = this.view.hGUI.tabOut;
            dataOut = get(tableOut, 'Data');
            
            if(col == IsotermTableRow.columnShow)
                hLine = this.view.axisHandles(row);
                hConf = this.view.hConfIntAxes(row, :);
                isChecked =  get(this.view.hGUI.cbConfInt, 'Value');
                
                if(dataOut{row, col} == true)
                    if ishandle(hLine)
                        set(hLine,'Visible','on');
                        if ishandle(hConf(1)) && ishandle(hConf(2))
                            set(hConf(1),'Visible','on');
                            set(hConf(2),'Visible','on');
                        end
                    end
                else
                    if ishandle(hLine)
                        set(hLine,'Visible','off');
                        if ishandle(hConf(1)) && ishandle(hConf(2)) 
                            set(hConf(1),'Visible','off');
                            set(hConf(2),'Visible','off');
                        end
                    end
                end
            end
        end
    end
    
    methods(Static = true)
        function copyExcel(o, e, table)
            str = clipboard('paste');
            str = strrep(str,',','.');
            import = str2num(str);
            set(table, 'Data', import);
        end
        
        function onClickPlusRow(o, e, table)
            data = get(table, 'Data');
            data(end+1,:) = [0, 0];
            set(table, 'Data', data);
        end
        
        function onClickMinusRow(o, e, table)
            data = get(table, 'Data');
            set(table, 'Data', data(1:end-1,:));
        end
        
    end
end
