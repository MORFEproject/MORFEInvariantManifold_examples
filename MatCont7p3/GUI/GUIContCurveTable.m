 classdef GUIContCurveTable
     
     
     
    methods(Static)
        function previewPanel(c, name, session)
            f = figure('Visible' , 'off' , 'NumberTitle' , 'off' , 'Name' , name, 'tag', 'matcont');
            if ~isempty(session)           
                listener = session.addlistener('shutdown', @(o, e) close(f));
            else
                listener = [];
            end
            mainnode = LayoutNode(-1,-1,'vertical');
            mainnode.setOptions('Add' , true);
            
            interpret = c.getStructInterpreter();
            
            h1 = GUIContCurveTable.makeTable(f, c.x , c.v , c.h , c.f, interpret);
            h2 = GUIContCurveTable.makeTable_s(f, c.s );
            
            s_tableresizer(h2);
            
            
            
            mainnode.addHandle(3 , 1 , h2 , 'minsize' , [Inf, 0] );
            
            h3 = uitable( f );
            mainnode.addHandle(3 , 1 , h3 , 'minsize' , [Inf, 100]); 
            
            set(h2 , 'CellSelectionCallback' , @(o,e) GUIContCurveTable.updateSdata(o,e,h3 , c.s) );
            mainnode.addHandle(10 , 1 , h1, 'minsize' , [Inf, Inf]);
            
            mainnode.makeLayoutHappen( get(f,'Position'));
            
            set(f,'Visible' , 'on' , 'ResizeFcn' , @(o,e)  GUIContCurveTable.layoutResize( o, e), 'DeleteFcn' , @(o,e) deleteobjects(mainnode, listener), 'UserData', timer('StartDelay', 0.2, 'TimerFcn', @(~, ~) mainnode.makeLayoutHappen( get(f,'Position')) ));
        end
        
        function layoutResize(o, e)
            stop(o.UserData);
            start(o.UserData);
            
        end
        
        
        function updateSdata(o,e,tablehandle , s_struct)
            GUIContCurveTable.fill_Sdata(tablehandle, s_struct( e.Indices(1)).data);
        end
        
        function handle = makeTable_s(parent,s,varargin)
            len = length(s);
            table = cell(len,2);
            for i = 1:len
                table{i,1} = s(i).label;
                table{i,2} = num2str(s(i).index);
                table{i,3} = s(i).msg;
            end
            handle = uitable(parent, 'Data' , table , 'RowName', [], 'ColumnName' , {'Label' , 'Index' , 'Message'}, varargin{:});
        end
        
        function handle = makeTable(parent, x,v,h,f, interpret, varargin)
            len = size(x,2);
            n = ones(1,len);
            blankline = mat2cell(blanks(len),1,n);
            
            m_x = size(x,1);
            table_x = mat2cell(x,ones(1,m_x) , n);
            
	    if(~isempty(v))
	    	m_v = size(v,1);
            	table_v = mat2cell(v,ones(1,m_v) , n);
            else
		m_v = 0;
		table_v = [];
	    end

	    if(~isempty(h))
		    m_h = size(h,1);
		    table_h = mat2cell(h,ones(1,m_h) , n);
            else
		    m_h = 0;
		    table_h = [];
	    end
        m_f = 0;
        table_f = [];
        try
            if (~isempty(f))
                m_f = size(f,1);
                table_f = mat2cell(f(:, 1:length(n)),ones(1,m_f) , n);
            end
        catch
        end
            rownames = {};
            for i=1:m_x
                rownames{end+1} = ['x(' num2str(i) ',:): ' interpret.x{i}];
            end
            rownames{end+1}  = '';
            for i=1:m_v
                rownames{end+1} = ['v(' num2str(i) ',:)'];
            end
            rownames{end+1}  = '';
            for i=1:m_h
                rownames{end+1} = ['h(' num2str(i) ',:): ' interpret.h{i}];
            end
            rownames{end+1}  = '';
            for i=1:m_f
                rownames{end+1} = ['f(' num2str(i) ',:): ' interpret.f{i}];
            end
            maxlen = max(cellfun(@(x) length(x), rownames));
            rownames = cellfun(@(x) {pad(x, maxlen)}, rownames);
            table = [table_x;blankline;table_v;blankline;table_h;blankline;table_f];
            handle = uitable(parent, 'Data' , table , 'RowName', rownames , varargin{:});
        end
        
        function fill_Sdata(tablehandle, sd)
           
            table = cell(0,0);
            rownames = cell(0,0);
            
            %FIXME: sd niet noodzakelijk altijd een struct , check
            
	    if (isstruct(sd))
	    fields = fieldnames(sd);
            
            
            row = 1;
            for k = 1:length(fields)
                rownames{row} = fields{k};
                val = sd.(fields{k});
                if (isempty( val))
                   table{row,1} = '[]';
                   row = row + 1;                    
                    
                    
                elseif (isnumeric( val ) )
                    [x,y] = size(val);
                    for i=1:x
                       for j = 1:y 
                         table{row + i - 1,j} = val(i,j);
                       end
                    end
                    row = row + x;
                    
                    
                else
                   table{row,1} = '???';
                   row = row + 1;
                end
            end
            
            
		    set(tablehandle, 'Data' , table , 'RowName', rownames  , 'ColumnName' , []);
            else
		    set(tablehandle, 'Data' , [] , 'RowName', []  , 'ColumnName' , []);
	    end
        end
        
        
        
        
        function test(s , i)
           openvar s(i).data; 
        end
    end
end

function s_tableresizer(handle)
    data = get(handle,'Data');
    col1 = maxChar(data(:,1)) * 7 + 40;
    col2 = maxChar(data(:,2)) * 7 + 40;
    col3 = maxChar(data(:,3)) * 7;
    
    set(handle , 'ColumnWidth' , {col1,col2,col3});


end
function m = maxChar(cellstr)
l = zeros(1, length(cellstr));
for i = 1:length(cellstr)
   l(i) = length(cellstr{i}); 
end
m = max(l);
end

function deleteobjects(o1, o2)
    delete(o1);
    delete(o2);

end
