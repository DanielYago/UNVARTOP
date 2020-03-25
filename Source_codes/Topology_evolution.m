%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
% A GUI to display topology evolution computed using any version of UNVARTOP code 	   %
%   Paper's Title                                                                      %
%		Topology Optimization using the UNsmooth VARiational Topology OPtimization 	   %
%			(UNVARTOP) method: an educational implementation in Matlab                 %
%   Authors                                                                            %
%       Daniel Yago, Juan Carlos Cante, Oriol Lloberas-Valls, Javier Oliver            %
% Versions                                                                             %
%   26/3/2020 - Submission version (Topology evolution)                                %
% GitHub repository																	   %
%	This code and other modifications of it can be downloaded from the Website:		   %
%		https://github.com/DanielYago/UNVARTOP	. Please send your comments to		   %
%		the author: daniel.yago@upc.edu								   				   %
% Copyright																			   %
%	Copyright (C) 2020  Daniel Yago													   %
%																					   %
%	This program is free software: you can redistribute it and/or modify it under the  %
%	terms of the GNU General Public License as published by the Free Software 		   %
%	Foundation, either version 3 of the License, or (at your option) any later version.%
%																					   %
%	This program is distributed in the hope that it will be useful, but WITHOUT ANY    %
%	WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    %
%	PARTICULAR PURPOSE.  See the GNU General Public License for more details.          %
%																					   %
%	You should have received a copy of the GNU General Public License along with this  %
%	program.  If not, see <https://www.gnu.org/licenses/>.							   %
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%%
classdef Topology_evolution < handle
    
    properties (Access=public)
        
        t_vec                           %vector of volumes
        coord                           %coordinates [nnodes,ndim]
        connect                         %connectivities [nel,nnodes_e]
        
        psi                             %psi nodal [nnodes,ntime]
        U                               %displacements [nnodes*nunkn,nloads,ntime]
        chi                             %chi nodal [nnodes,ntime]
        
        coord_mapping                   %coordinates for the corresponding t_ref
        connect_mapping                 %connectivities for the corresponding t_ref
        U_mapping                       %displacements for the corresponding t_ref
        
    end
    properties (Access=protected)
        
        chi_e                           %chi elemental [nel,ntime]
        
        t_ref               = 1;        %time reference to plot
        scale_elem          = 10;       %scale of displacements (normalized with maximum number of elements)
        scale               = 1;        %scale of displacements (deformed mesh)
        dt                  = 0.1;      %dt between steps (for animation)
        
        fg_handle                       %figure handle
        
        ax1_handle                      %handle of axes 1
        ax2_handle                      %handle of axes 2
        
        obj_ax1_handle                  %handle of object in axes 1
        obj_ax2_fix_handle              %handle of object 1 in axes 2
        obj_ax2_free_handle             %handle of object 2 in axes 2
        
        pbt_play                        %pushbutton play
        pbt_pause                       %pushbutton pause
        pbt_backward_1                  %pushbutton move backward 1 step
        pbt_backward_end                %pushbutton go to first
        pbt_forward_1                   %pushbutton move forward 1 step
        pbt_forward_end                 %pushbutton go to last
        tgbt_loop                       %pushbutton loop
        txt_loop                        %text for loop
        dd_looptype                     %drop-down for loop type
        
        txt_color                       %text for color
        dd_color                        %drop-down for color selection
        txt_style                       %text for style
        dd_style                        %drop-down for style selection
        txt_fscale                      %text for field to scale
        dd_fscale                       %drop-down for field to scale selection
        txt_scale                       %text for scale
        txtedit_scale                   %text area for scale input
        txt_time                        %text for time between steps
        txtedit_time                    %text are for time input
        
        txt_sym                         %text for symmetries
        chb_xmin                        %checkbox xmin symmetry
        chb_xmax                        %checkbox xmax symmetry
        chb_ymin                        %checkbox ymin symmetry
        chb_ymax                        %checkbox ymax symmetry
        
        sym_x_min                       %boolean for xmin symmetry
        sym_x_max                       %boolean for xmax symmetry
        sym_y_min                       %boolean for ymin symmetry
        sym_y_max                       %boolean for ymax symmetry
        
        is_animated         = false;    %boolean to know if it is animated
        
        %faces_mapping       = [1 2 5;2 3 5;3 4 5;4 1 5]'; %mesh mapping for plotting nodal fields (quadrilateral mesh -> 4*triangular mesh)
        faces_mapping       = [1 2 3 4]'; %mesh mapping for plotting nodal fields
    end
    
    methods
        function obj = Topology_evolution(varargin)
            
            obj.coord = varargin{1};
            obj.connect = varargin{2};
            obj.t_vec = varargin{3};
            
            obj.psi = varargin{4};
            obj.chi_e = varargin{5};
            obj.U = varargin{6};
            
            obj.smooth_chi;
            
            obj.compute_variable_mapping;
            obj.detect_symmetries;
            
            obj.create_figure;
        end
        
        function smooth_chi(obj)
            % function that smooth the elemental chi field using a Global
            % smoothing (Mass matrix)
            
            ME_Lap = 1/36*[ 4  2  1  2; 2 4  2  1; 1  2 4  2; 2  1  2 4];
            i_KF = reshape(kron(obj.connect,ones(4,1))',4^2*size(obj.connect,1),1);
            j_KF = reshape(kron(obj.connect,ones(1,4))',4^2*size(obj.connect,1),1);
            s_KF = reshape(ME_Lap(:)*ones(1,size(obj.connect,1)),4^2*size(obj.connect,1),1);
            M_Lap = sparse(i_KF,j_KF,s_KF);
            
            posgp1 = [0;0];
            W1 = 4;
            i_xi = reshape(obj.connect',4*size(obj.connect,1),1);
            N_matrix = 0.25*(1+[-1 1 1 -1]'*posgp1(1,:)).*(1+[-1 -1 1 1]'*posgp1(2,:));
            N_T = N_matrix.*W1/4;
            
            obj.chi=M_Lap\accumarray([kron(ones(numel(obj.t_vec),1),i_xi) kron((1:numel(obj.t_vec))',ones(size(i_xi)))],...
                reshape(kron(obj.chi_e(:,1:numel(obj.t_vec)),N_T),[],1),[size(obj.coord,1) numel(obj.t_vec)]);
            
        end
        
        function create_figure(obj)
            
            %% Create figure
            
            obj.fg_handle = figure();
            set(obj.fg_handle,'Name','Topology Evolution','Color','w','Units','normalized','Resize','off');
            
            %% Create axes
            
            obj.ax1_handle = axes(obj.fg_handle,'Position',[0.1 0.125 0.8 0.625],'Box','off','XTick', [],'YTick', []);
            axis(obj.ax1_handle,'equal','tight');
            obj.ax2_handle = axes(obj.fg_handle,'Position',[0.1 0.05 0.8 0.05],'Box','on','YTick', []);
            
            %% Create buttons sequence
            
            obj.pbt_play = uicontrol(obj.fg_handle,'Style','pushbutton','String',sprintf('\x25b6'),...
                'Units','normalized','Position',[0.1 0.875 0.06 0.06],...
                'Callback',@(src,event)obj.animate_topology);
            
            obj.pbt_pause = uicontrol(obj.fg_handle,'Style','pushbutton','String','| |',...
                'Units','normalized','Position',[0.175 0.875 0.06 0.06],...
                'Callback','uiwait(gcbf)');
            
            obj.pbt_backward_end = uicontrol(obj.fg_handle,'Style','pushbutton','String','|<',...
                'Units','normalized','Position',[0.25 0.875 0.06 0.06],...
                'Callback',@(src,event)obj.advance_time(-numel(obj.t_vec)));
            
            obj.pbt_backward_1 = uicontrol(obj.fg_handle,'Style','pushbutton','String','<',...
                'Units','normalized','Position',[0.325 0.875 0.06 0.06],...
                'Callback',@(src,event)obj.advance_time(-1));
            
            obj.pbt_forward_1 = uicontrol(obj.fg_handle,'Style','pushbutton','String','>',...
                'Units','normalized','Position',[0.4 0.875 0.06 0.06],...
                'Callback',@(src,event)obj.advance_time(1));
            
            obj.pbt_forward_end = uicontrol(obj.fg_handle,'Style','pushbutton','String','>|',...
                'Units','normalized','Position',[0.475 0.875 0.06 0.06],...
                'Callback',@(src,event)obj.advance_time(numel(obj.t_vec)));
            
            obj.tgbt_loop = uicontrol(obj.fg_handle,'Style','togglebutton','String',sprintf('\x27f3'),...
                'Units','normalized','Position',[0.55 0.875 0.06 0.06],...
                'Callback',@(src,event)obj.animate_loop_topology);
            
            obj.txt_loop = uicontrol(obj.fg_handle,'Style','text',...
                'Units','normalized','Position',[0.625 0.855 0.1 0.06],...
                'String','Loop stl:','HorizontalAlignment','left',...
                'BackgroundColor','w');
            
            obj.dd_looptype = uicontrol(obj.fg_handle,'Style','popupmenu',...
                'Units','normalized','Position',[0.735 0.855 0.165 0.075],...
                'String',{'Volume','Scale linear','Scale sine'},'Value',1);
            
            %% Create buttons Color, style, scale and time
            
            obj.txt_color = uicontrol(obj.fg_handle,'Style','text',...
                'Units','normalized','Position',[0.1 0.79 0.1 0.06],...
                'String','Color:','HorizontalAlignment','left',...
                'BackgroundColor','w');
            
            temp = strcat('U',num2str((1:size(obj.U,2))')); temp = cellstr(temp)';
            obj.dd_color = uicontrol(obj.fg_handle,'Style','popupmenu',...
                'Units','normalized','Position',[0.165 0.79 0.07 0.075],...
                'String',[{'Psi'},{'Chi'},temp(:)'],'Value',1,...
                'Callback',@(src,event)obj.change_color);
            
            obj.txt_style = uicontrol(obj.fg_handle,'Style','text',...
                'Units','normalized','Position',[0.25 0.79 0.1 0.06],...
                'String','Style:','HorizontalAlignment','left',...
                'BackgroundColor','w');
            
            obj.dd_style = uicontrol(obj.fg_handle,'Style','popupmenu',...
                'Units','normalized','Position',[0.315 0.79 0.08 0.075],...
                'String',{'Surf','Wire','Wire&Surf'},'Value',1,...
                'Callback',@(src,event)obj.update_style);
            
            obj.txt_fscale = uicontrol(obj.fg_handle,'Style','text',...
                'Units','normalized','Position',[0.4 0.79 0.1 0.06],...
                'String','Displ:','HorizontalAlignment','left',...
                'BackgroundColor','w');
            
            obj.dd_fscale = uicontrol(obj.fg_handle,'Style','popupmenu',...
                'Units','normalized','Position',[0.46 0.79 0.075 0.075],...
                'String',temp(:)','Value',1,...
                'Callback',@(src,event)obj.update_field_displacement);
            
            obj.txt_scale = uicontrol(obj.fg_handle,'Style','text',...
                'Units','normalized','Position',[0.55 0.79 0.175 0.06],...
                'String','Scale[el]:','HorizontalAlignment','left',...
                'BackgroundColor','w');
            
            obj.txtedit_scale = uicontrol(obj.fg_handle,'Style','edit',...
                'Units','normalized','Position',[0.64 0.8 0.075 0.06],...
                'String','1',...
                'Callback',@(src,event)obj.update_scale);
            
            obj.txt_time = uicontrol(obj.fg_handle,'Style','text',...
                'Units','normalized','Position',[0.72 0.79 0.15 0.06],...
                'String','Time[s]:','HorizontalAlignment','left',...
                'BackgroundColor','w');
            
            obj.txtedit_time = uicontrol(obj.fg_handle,'Style','edit',...
                'Units','normalized','Position',[0.8 0.8 0.1 0.06],...
                'String','0.1',...
                'Callback',@(src,event)obj.update_time);
            
            %% Symmetry buttons
            
            obj.txt_sym = uicontrol(obj.fg_handle,'Style','text',...
                'Units','normalized','Position',[0.1 0.74 0.1 0.06],...
                'String','Sym:','HorizontalAlignment','left',...
                'BackgroundColor','w');
            
            obj.chb_xmin = uicontrol(obj.fg_handle,'Style','checkbox','String','xmin',...
                'Units','normalized','Position',[0.165 0.75 0.1 0.06],...
                'BackgroundColor','w','Visible',obj.sym_x_min,...
                'Callback',@(src,event)obj.create_symmetry_xmin);
            
            obj.chb_xmax = uicontrol(obj.fg_handle,'Style','checkbox','String','xmax',...
                'Units','normalized','Position',[0.275 0.75 0.1 0.06],...
                'BackgroundColor','w','Visible',obj.sym_x_max,...
                'Callback',@(src,event)obj.create_symmetry_xmax);
            
            obj.chb_ymin = uicontrol(obj.fg_handle,'Style','checkbox','String','ymin',...
                'Units','normalized','Position',[0.385 0.75 0.1 0.06],...
                'BackgroundColor','w','Visible',obj.sym_y_min,...
                'Callback',@(src,event)obj.create_symmetry_ymin);
            
            obj.chb_ymax = uicontrol(obj.fg_handle,'Style','checkbox','String','ymax',...
                'Units','normalized','Position',[0.495 0.75 0.1 0.06],...
                'BackgroundColor','w','Visible',obj.sym_y_max,...
                'Callback',@(src,event)obj.create_symmetry_ymax);
            
            %% Print Topology and volume time-line
            
            obj.obj_ax1_handle = patch(obj.ax1_handle,'Vertices',obj.coord_mapping,'Faces',obj.connect_mapping,'EdgeColor','none','FaceColor','interp');
            obj.update_scale;
            obj.change_color;
            
            obj.obj_ax2_fix_handle = patch(obj.ax2_handle,'Vertices',[repelem(obj.t_vec',2,1), repmat([0; 1],numel(obj.t_vec),1)],...
                'Faces',reshape(1:2*numel(obj.t_vec),2,[])','EdgeColor','blue','FaceColor','blue');
            obj.obj_ax2_free_handle = patch(obj.ax2_handle,'Vertices',[repmat(obj.t_vec(obj.t_ref),2,1), [0; 1]],'Faces',[1 2],'EdgeColor','red','FaceColor','red','LineWidth',4);
            
        end
        
        function animate_topology(obj)
            % Callback function when pbt_play is pushed. It starts the
            % animation in the current t_ref till the last one, modifying the
            % volume. It can be paused by pushing the pbt_pause button.
            
            if ~obj.is_animated
                while obj.t_ref < numel(obj.t_vec)
                    obj.is_animated = true;
                    advance_time(obj,1);
                    drawnow;
                    pause(obj.dt);
                end
                obj.is_animated = false;
            else
                uiresume(obj.fg_handle);
                obj.is_animated = false;
            end
        end
        
        function animate_loop_topology(obj)
            % Callback function when tgbt_loop is pushed. Depending on the
            % value of obj.dd_looptype, the loop is done throughout the volume
            % ratio (t_vec) or for a given t_vec, the scale factor of the
            % displacements is modified, in a linear o sine way.
            
            scl = 0;
            if get(obj.tgbt_loop,'Value')
                while get(obj.tgbt_loop,'Value')
                    switch get(obj.dd_looptype,'Value')
                        case 1
                            advance_time(obj,1);
                            if obj.t_ref == numel(obj.t_vec)
                                obj.t_ref = 0;
                            end
                        case 2
                            scl = scl + 1/50;
                            advance_scale(obj,scl);
                            if abs(scl-1)<1e-6
                                scl = 0;
                            end
                        case 3
                            scl = scl + 1/50;
                            advance_scale(obj,sin(2*pi*scl));
                            if abs(scl-1)<1e-6
                                scl = 0;
                            end
                    end
                    drawnow;
                    pause(obj.dt);
                end
            end
            
        end
        
        function advance_time(obj,step)
            % function that increases t_ref in step and shows the corresponding plot
            % for the new t_ref for obj_ax1_handle and obj_ax2_free_handle
            
            obj.t_ref = max(min(obj.t_ref+step,numel(obj.t_vec)),1);
            
            obj.update_variable_mapping;
            obj.apply_symmetries;
            obj.change_color;
            obj.update_coord(obj.coord_mapping+obj.scale*obj.U_mapping);
            
            set(obj.obj_ax2_free_handle,'Vertices',[repmat(obj.t_vec(obj.t_ref),2,1), [0; 1]]);
            
        end
        
        function advance_scale(obj,scl_factor)
            % function that modifies the scale of the displacements. It is used
            % for the scale-displacement loops
            
            fdispl = get(obj.dd_fscale,'Value'); if isempty(fdispl); fdispl=1; end
            scl_max = obj.scale_elem/max(abs(reshape(obj.U(:,fdispl,obj.t_ref),[],1)));
            
            obj.update_coord(obj.coord_mapping+scl_factor*scl_max*obj.U_mapping);
        end
        
        function change_color(obj)
            % Callback function when dd_color is modified. It updates the color
            % of the plot for the corresponding mesh mapping and t_ref
            
            switch get(obj.dd_color,'Value')
                case 1 %'Psi'
                    Color_mapping = [obj.psi(:,obj.t_ref); mean(reshape(obj.psi(obj.connect',obj.t_ref),4,[]),1)'];
                    caxis(obj.ax1_handle,[-1 1]); colormap(obj.ax1_handle,flip(gray(2)));
                case 2 %'Chi'
                    Color_mapping = [obj.chi(:,obj.t_ref); mean(reshape(obj.chi(obj.connect',obj.t_ref),4,[]),1)'];
                    caxis(obj.ax1_handle,[0 1]); colormap(obj.ax1_handle,flip(gray(20)));
                otherwise %'U'
                    fdispl = get(obj.dd_color,'Value')-2;
                    u_rshp = reshape(obj.U(:,fdispl,obj.t_ref),2,[])'; u_norm = sqrt(sum(u_rshp.^2,2));
                    Color_mapping = [u_norm; mean(reshape(u_norm(obj.connect'),4,[]),1)'];
                    caxis(obj.ax1_handle,[min(Color_mapping(:)) max(Color_mapping(:))]); colormap(obj.ax1_handle,jet(20));
            end
            
            n_sym = get(obj.chb_xmin,'Value')+get(obj.chb_xmax,'Value')+get(obj.chb_ymin,'Value')+get(obj.chb_ymax,'Value');
            Color_mapping = [Color_mapping;repmat(Color_mapping,n_sym,1)];
            
            obj.update_color(Color_mapping);
            
        end
        
        function compute_variable_mapping(obj)
            % function that updates the mapping coordinates, connectivities and
            % displacement for the corresponding t_ref (t_ref for the initial
            % iteration is equal to 1)
            
            update_variable_mapping(obj);
            
            obj.coord_mapping = [obj.coord; mean(reshape(obj.coord(obj.connect',1),4,[]),1)' mean(reshape(obj.coord(obj.connect',2),4,[]),1)'];
            
            obj.connect_mapping = [obj.connect size(obj.coord,1)+(1:size(obj.connect,1))'];
            obj.connect_mapping = reshape(obj.connect_mapping(:,obj.faces_mapping(:))',size(obj.faces_mapping,1),[])';
        end
        
        function update_variable_mapping(obj)
            % function that updates the mapping displacements for the
            % corresponding t_ref
            
            fdispl = get(obj.dd_fscale,'Value'); if isempty(fdispl); fdispl=1; end
            u_rshp = reshape(obj.U(:,fdispl,obj.t_ref),2,[])';
            obj.U_mapping = [u_rshp; mean(reshape(u_rshp(obj.connect',1),4,[]),1)' mean(reshape(u_rshp(obj.connect',2),4,[]),1)'];
            
        end
        
        function update_coord(obj,coord_new)
            % function that updates the coordinates of obj_ax1_handle
            
            set(obj.obj_ax1_handle,'Vertices',coord_new);
            
        end
        
        function update_connect(obj,connect_new)
            % function that updates the connectivities of obj_ax1_handle
            
            set(obj.obj_ax1_handle,'Faces',connect_new);
            
        end
        
        function update_color(obj,color_new)
            % function that updates the color of obj_ax1_handle
            
            set(obj.obj_ax1_handle,'FaceVertexCData',color_new);
            
        end
        
        function update_scale (obj)
            % Callback function when txtedit_scale is edited. It updates the scale factor
            % of the displacements (modifying the coordinates of obj_ax1_handle)
            
            scl = str2double(get(obj.txtedit_scale,'String'));
            
            if ~ischar(abs(scl))
                obj.scale_elem = scl;
                fdispl = get(obj.dd_fscale,'Value'); if isempty(fdispl); fdispl=1; end
                obj.scale = obj.scale_elem/max(abs(reshape(obj.U(:,fdispl,1:end),[],1)));
                
                obj.update_coord(obj.coord_mapping+obj.scale*obj.U_mapping);
            else
                set(obj.txtedit_scale,'String',num2str(obj.scale_elem));
            end
        end
        
        function update_style(obj)
            % Callback function when dd_style is modified. It updates the style or
            % properties of the patch (obj_ax1_handle)
            
            switch get(obj.dd_style,'Value')
                case 1 %Surface
                    set(obj.obj_ax1_handle,'EdgeColor','none','FaceColor','interp');
                case 2 %Wireframe
                    set(obj.obj_ax1_handle,'EdgeColor','blue','FaceColor','none');
                case 3 %Wireframe and surface
                    set(obj.obj_ax1_handle,'EdgeColor','blue','FaceColor','interp');
                    obj.change_color;
            end
            
        end
        
        function update_field_displacement(obj)
            % Callback function when dd_fscale is modified. It updates the
            % displacement field which is used to deform the mesh of the patch (obj_ax1_handle)
            
            obj.update_variable_mapping;
            obj.update_scale;
            
        end
        
        function create_symmetry_xmin(obj)
            % Callback function when checkbox chb_xmin is checked. It mirrors
            % coordiantes, connectivities and displacement with respect
            % to the minimum x coordinate
            
            if get(obj.chb_xmin,'Value')
                n_coord = size(obj.coord_mapping,1);
                obj.coord_mapping = [obj.coord_mapping; 2*min(obj.coord_mapping(:,1))-obj.coord_mapping(:,1)  obj.coord_mapping(:,2)];
                obj.connect_mapping = [obj.connect_mapping; n_coord+obj.connect_mapping];
                obj.U_mapping = [obj.U_mapping; -obj.U_mapping(:,1) obj.U_mapping(:,2)];
                
                obj.update_connect(obj.connect_mapping);
                obj.update_scale;
                obj.change_color;
            else
                obj.compute_variable_mapping;
                
                obj.update_connect(obj.connect_mapping);
                obj.update_scale;
                obj.change_color;
            end
            
        end
        
        function create_symmetry_xmax(obj)
            % Callback function when checkbox chb_xmax is checked. It  mirrors
            % coordiantes, connectivities and displacement with respect
            % to the maximum x coordinate
            
            if get(obj.chb_xmax,'Value')
                n_coord = size(obj.coord_mapping,1);
                obj.coord_mapping = [obj.coord_mapping; 2*max(obj.coord_mapping(:,1))-obj.coord_mapping(:,1)  obj.coord_mapping(:,2)];
                obj.connect_mapping = [obj.connect_mapping; n_coord+obj.connect_mapping];
                obj.U_mapping = [obj.U_mapping; -obj.U_mapping(:,1) obj.U_mapping(:,2)];
                
                obj.update_connect(obj.connect_mapping);
                obj.update_scale;
                obj.change_color;
            else
                obj.compute_variable_mapping;
                
                obj.update_connect(obj.connect_mapping);
                obj.update_scale;
                obj.change_color;
            end
            
        end
        
        function create_symmetry_ymin(obj)
            % Callback function when checkbox chb_ymin is checked. It mirrors
            % coordiantes, connectivities and displacement with respect
            % to the minimum y coordinate
            
            if get(obj.chb_ymin,'Value')
                n_coord = size(obj.coord_mapping,1);
                obj.coord_mapping = [obj.coord_mapping; obj.coord_mapping(:,1) 2*min(obj.coord_mapping(:,2))-obj.coord_mapping(:,2) ];
                obj.connect_mapping = [obj.connect_mapping; n_coord+obj.connect_mapping];
                obj.U_mapping = [obj.U_mapping; obj.U_mapping(:,1) -obj.U_mapping(:,2)];
                
                obj.update_connect(obj.connect_mapping);
                obj.update_scale;
                obj.change_color;
            else
                obj.compute_variable_mapping;
                
                obj.update_connect(obj.connect_mapping);
                obj.update_scale;
                obj.change_color;
            end
            
        end
        
        function create_symmetry_ymax(obj)
            % Callback function when checkbox chb_ymax is checked. It mirrors
            % coordiantes, connectivities and displacement with respect
            % to the maximum y coordinate
            
            if get(obj.chb_ymax,'Value')
                n_coord = size(obj.coord_mapping,1);
                obj.coord_mapping = [obj.coord_mapping; obj.coord_mapping(:,1) 2*max(obj.coord_mapping(:,2))-obj.coord_mapping(:,2)];
                obj.connect_mapping = [obj.connect_mapping; n_coord+obj.connect_mapping];
                obj.U_mapping = [obj.U_mapping; obj.U_mapping(:,1) -obj.U_mapping(:,2)];
                
                obj.update_connect(obj.connect_mapping);
                obj.update_scale;
                obj.change_color;
            else
                obj.compute_variable_mapping;
                
                obj.update_connect(obj.connect_mapping);
                obj.update_scale;
                obj.change_color;
            end
            
        end
        
        function apply_symmetry_xmin(obj)
            % function that applies the symmetry on the displacement mapping
            % when the symmetry with respect to the minimum x coordinates is
            % applied
            
            if get(obj.chb_xmin,'Value')
                obj.U_mapping = [obj.U_mapping; -obj.U_mapping(:,1) obj.U_mapping(:,2)];
            end
            
        end
        
        function apply_symmetry_xmax(obj)
            % function that applies the symmetry on the displacement mapping
            % when the symmetry with respect to the maximum x coordinates is
            % applied
            
            if get(obj.chb_xmax,'Value')
                obj.U_mapping = [obj.U_mapping; -obj.U_mapping(:,1) obj.U_mapping(:,2)];
            end
            
        end
        
        function apply_symmetry_ymin(obj)
            % function that applies the symmetry on the displacement mapping
            % when the symmetry with respect to the minimum y coordinates is
            % applied
            
            if get(obj.chb_ymin,'Value')
                obj.U_mapping = [obj.U_mapping; obj.U_mapping(:,1) -obj.U_mapping(:,2)];
            end
            
        end
        
        function apply_symmetry_ymax(obj)
            % function that applies the symmetry on the displacement mapping
            % when the symmetry with respect to the maximum y coordinates is
            % applied
            
            if get(obj.chb_ymax,'Value')
                obj.U_mapping = [obj.U_mapping; obj.U_mapping(:,1) -obj.U_mapping(:,2)];
            end
            
        end
        
        function apply_symmetries(obj)
            
            obj.apply_symmetry_xmin;
            obj.apply_symmetry_xmax;
            obj.apply_symmetry_ymin;
            obj.apply_symmetry_ymax;
            
        end
        
        function detect_symmetries(obj)
            % function that detects which of the boundary sides are possible
            % symmetry axis
            
            id_x_min = obj.coord(:,1)==min(obj.coord(:,1));
            id_x_max = obj.coord(:,1)==max(obj.coord(:,1));
            id_y_min = obj.coord(:,2)==min(obj.coord(:,2));
            id_y_max = obj.coord(:,2)==max(obj.coord(:,2));
            
            u_rshp = reshape(obj.U(:,1,1),2,[])';
            
            obj.sym_x_min = ~any(u_rshp(id_x_min,1));
            obj.sym_x_max = ~any(u_rshp(id_x_max,1));
            obj.sym_y_min = ~any(u_rshp(id_y_min,2));
            obj.sym_y_max = ~any(u_rshp(id_y_max,2));
            
        end
        
        function update_time(obj)
            % Callback function when txtedit_time is edited. It modifies the dt
            % between steps, when it is animated
            
            dtime = str2double(get(obj.txtedit_time,'String'));
            
            if ~ischar(abs(dtime))
                obj.dt = abs(dtime);
            else
                set(obj.txtedit_time,'String',num2str(obj.dt));
            end
            
        end
        
    end
    
end