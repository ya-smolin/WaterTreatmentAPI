classdef IsotermTableRow
    properties(Constant)
        columnIsotermName = 1;
        columnShow = 4;
        columnConfInt = 7;
        columnName =      {'name',    'formula',  'koef',   'show',     'rsquare', 'sse', 'koef range'};
        columnFormat =    {'char',    'char',     'char',   'logical',  'numeric', 'numeric', 'char'}
        columnEditable =  [false,   false,    false,  true,     false,   false,   false]
        size = length(IsotermTableRow.columnName);
        CONFIDENCE_LEVEL_DEFAULT = 0.95;
    end
    properties
        data;
        axisHandles = -1;
        hConfIntAxes = -ones(1, 2);
    end
    
    methods
        function initHoverLegend(this)
            if(ishghandle(this.axisHandles))
                hMenu = uicontextmenu;
                hMenuItem = uimenu(hMenu, 'Label', this.data{IsotermTableRow.columnIsotermName}, 'HandleVisibility','callback');
                set(this.axisHandles, 'uicontextmenu',hMenu);
                 
            end
        end
        
        function updatePlot(this, isConfidencePlotTurnOn)
             if(this.isPlotVisible())
                    this.visibleAxes(true);
                    if(isConfidencePlotTurnOn)
                        this.visibleAxesConf(true);
                    else
                        this.visibleAxesConf(false);
                    end
                else
                    this.visibleAxes(false);
                    this.visibleAxesConf(false);
             end
        end
        
        function isVisible = isPlotVisible(this)
            isVisible = this.data{IsotermTableRow.columnShow};
        end
        
        function visibleAxes(this, turnOn)
            if this.axisHandles ~= -1 && ishghandle(this.axisHandles)
                    if(turnOn)
                        set(this.axisHandles,'Visible','on');
                    else
                        set(this.axisHandles,'Visible','off');
                    end
            end
        end
        
        function clearPlots(this)
            plotLine = this.axisHandles;
            if(plotLine ~= -1 && ishghandle(plotLine))
                    delete(plotLine);
            end;
            plotLineConf = this.hConfIntAxes;
            if sum(plotLineConf ~= -1 & ishghandle(plotLineConf)) == 2
                    delete(plotLineConf);
            end;
        end
        
        function visibleAxesConf(this, turnOn)
            if sum(this.hConfIntAxes ~= -1 & ishghandle(this.hConfIntAxes)) == 2
                if(turnOn)
                    set(this.hConfIntAxes,'Visible','on');
                else
                    set(this.hConfIntAxes,'Visible','off');
                end
            end
        end
            
        function this = IsotermTableRow(isoterm, confidenceLevel)
            if nargin == 1
                confidenceLevel = IsotermTableRow.CONFIDENCE_LEVEL_DEFAULT;
            end
            if isa(isoterm, 'IsotermType') 
                this.data = {
                    isoterm.name,...
                    isoterm.formula,...
                    '',...
                    false,...
                    0,...
                    Inf,...
                    ''
                    };
            else
                isotermType = isoterm.isotermType;
                isotermResult =  isoterm.isotermResult;
                gof = isoterm.gof;
                confInt = confint(isotermResult, confidenceLevel)';
                confInt(confInt<0) = 0;
                this.data = {
                    isotermType.name,...
                    isotermType.formula,...
                    mat2str(round(coeffvalues(isotermResult) * 100) / 100),...
                    true,...
                    gof.adjrsquare,...
                    gof.sse,...
                    mat2str(round(confInt * 100) / 100)
                    };
            end
        end
    end

end

