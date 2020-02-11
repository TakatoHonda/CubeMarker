function [result_f] = plots_results( dir, varargin )
    %% maximum # of sequences
    seq_max = 4;
    need_legend = false;
    plot_seg = false;

    %% load original data and segmenst information
    fid = fopen([dir, 'input']);
    fn = textscan(fid, '%s');
    fclose(fid);
    seq_length = length(fn{1});
    if seq_length > seq_max 
        seq_length = seq_max;
    end
    seq_list = (0:seq_length-1);
    
    %% load segment information
    datfn=[dir, 'segment.'];
    fid = fopen([datfn, 'labels']);
    labels = textscan(fid, '%s %s %s %s %s');
    n_regime = size(labels{1},1);  % number of regimes
    fclose(fid);
    
    %% main figure
    result_f = figure('Position', [0 0 800 400]);
    
    %% plot result for each sequence
    regime_list = zeros(1, n_regime);
    
    %% load original data
    orgfn=fn{1}{1};
    org = load(orgfn);
    for i_seq = 2:seq_length
        seq_id = seq_list(i_seq);
        %% select single sequence
        orgfn=fn{1}{seq_id+1};

        %% load orginal data
        org(:,:,i_seq) = load(orgfn);
    end
    
    %% normalization
    [n,d,w] = size(org);
    hoge = permute(org, [2,1,3]);  % d n w
    hoge = reshape(hoge, d, [])';
    hoge = zscore(hoge);  % d * nw
    org = reshape(hoge', d,n,w); 
    org = permute(org, [2,1,3]);
  
    st=1;
    ed=length(org)-1;
    miny=min(min(min(org)));
    maxy=max(max(max(org)));
    
    % plot each sequence
    for i_seq = 1:seq_length
        seq_id = seq_list(i_seq);
        figure(result_f);

        %% plot original data
        if seq_length == i_seq
            ax = subplot(seq_length+1, 1, i_seq, 'Position', [ 0.4-i_seq*0.2/seq_length 0.9-i_seq*0.5/seq_length 0.5 1/seq_length]);
        else
            ax = subplot(seq_length+1, 1, i_seq, 'Position', [ 0.4-i_seq*0.2/seq_length 0.9-i_seq*0.5/seq_length 0.5 1/seq_length]);
        end

        %% plot segments
        if plot_seg
            new_regime_list = plot_segments(datfn, n_regime, st, ed, miny, maxy, seq_id, seq_length, i_seq, org(:,:,i_seq), regime_list);
            regime_list = new_regime_list;
        end

        %% plot single sequence
        plot(org(:,:,i_seq), 'LineWidth', 1.0);
        xlim([st, ed])
        ylim([miny*1.1 maxy*1.1])
        yticklabels([]);
    end

    if need_legend
        %% make legend title
        title_list{1, n_regime} = [];
        for i = 1:n_regime
            mnum = num2str(i);
            title_list{1,i} = strcat('\theta', mnum);
        end
        %% show legend
        ax = subplot(seq_length+1, 1, seq_length+1);
        set(ax, 'Visible', 'off');
        legend(ax, regime_list, 'DisplayName', title_list, 'Orientation', 'horizontal', 'Location', 'southeast');
    end
    
    %% plot on google map
    if nargin ~= 1
        map_data = load(fullfile(varargin{1}));
        tmp.d = mat2dataset(map_data.X,'VarNames',map_data.md.col_names);
        map_f = plot_gmap(datfn, n_regime, tmp.d, seq_length);
    end
        
end


function [result_f map_f] = plots_each(dir, varargin)
    %% maximum # of sequences
    seq_max = 32;

    %% load original data and segmenst information
    fid = fopen([dir, 'input']);
    fn = textscan(fid, '%s');
    fclose(fid);
    seq_length = length(fn{1});
    if seq_length > seq_max 
        seq_length = seq_max;
    end
    seq_list = (0:seq_length-1);
    
    %% load segment information
    datfn=[dir, 'segment.'];
    fid = fopen([datfn, 'labels']);
    labels = textscan(fid, '%s %s %s %s %s');
    n_regime = size(labels{1},1);  % number of regimes
    fclose(fid);
    
    %% main figure
    result_f = figure('Position', [0 0 800 400]);
    
    %% plot result for each sequence
    regime_list = zeros(1, n_regime);
    
    %% load original data
    orgfn=fn{1}{1};
    org = load(orgfn);
    org_unfold = load(orgfn);
    for i_seq = 2:seq_length
        seq_id = seq_list(i_seq);
        %% select single sequence
        orgfn=fn{1}{seq_id+1};

        %% load orginal data
        org(:,:,i_seq) = load(orgfn);
        org_unfold = [org_unfold;load(orgfn)];
    end
    
    %% normalization
    [n,d,w] = size(org);
    hoge = permute(org, [2,1,3]);  % d n w
    hoge = reshape(hoge, d, [])';
    hoge = zscore(hoge);  % d * nw
    org = reshape(hoge', d,n,w); 
    org = permute(org, [2,1,3]);
  
    st=1;
    ed=length(org)-1;
%     miny=min(min(min(org(st:500))));
%     maxy=max(max(max(org(st:500))));
    miny=min(min(min(org)));
    maxy=max(max(max(org)));
    
    for i_seq = 1:seq_length
        seq_id = seq_list(i_seq);
        %% plot segments
        figure(result_f);
        [new_regime_list, ax] = plot_segments(datfn, n_regime, st, ed, miny, maxy, seq_id, seq_length, i_seq, org(:,:,i_seq), regime_list);
        regime_list = new_regime_list;
    end
    
    %% make legend title
    title_list{1, n_regime} = [];
    for i = 1:n_regime
        mnum = num2str(i);
        title_list{1,i} = strcat('\theta', mnum);
    end
    
    %% show legend
    ax = subplot(seq_length+1, 1, seq_length+1);
    set(ax, 'Visible', 'off');
    legend(ax, regime_list, 'DisplayName', title_list, 'Orientation', 'horizontal', 'Location', 'southeast');
        
    %% plot on google map
    if nargin ~= 1
        map_data = load(fullfile(varargin{1}));
        tmp.d = mat2dataset(map_data.X,'VarNames',map_data.md.col_names);
        map_f = plot_gmap(datfn, n_regime, tmp.d, seq_length);
    end
end


function mapf = plot_gmap(datfn, n_regime, d, seq_length)
    
    %% regime color and transparency setting
    alpha = 0.6*1/seq_length;
%     alpha = 0.05
    max_colors = 8;
    if max_colors > n_regime
        max_colors = n_regime;
    end

    orgf=gcf;
    mapf=figure;
    xlim([min(d.lon) max(d.lon)]);
    ylim([min(d.lat) max(d.lat)]);
    plot_google_map('language','ja','MapType','satellite');
    
    for i=0: n_regime-1
        figure(orgf);
        
        %% regime color setting 
        cols = colormap(jet(max_colors));
        col=cols(mod(i, size(cols, 1)) + 1,:);
%         if mod(i, 2) == 0
%             col = cols(i/2 + 1,:);
%         else
%             col = cols(size(cols, 1) - ceil(i/2) + 1, :);
%         end
        col = [col, alpha];
        
        %% load data: "segment.[i]"
        fid = fopen([datfn, num2str(i)]);
        tline = fgets(fid);
        X = [ ];
        while ischar(tline)
            dline = str2num(tline);
            % if dline(3) == seq_id
            if 1
                X = [X; dline];
            end
            tline = fgets(fid);
        end
        fclose(fid);
        
        figure(mapf);
        for j=1:size(X,1)
            p = plot(d.lon(X(j,1)+1:X(j,2)+1), d.lat(X(j,1)+1:X(j,2)+1), 'Color', col, 'LineWidth',5);
        end
    end
end


function [new_regime_list, ax] = plot_segments(datfn, n_regime, st, ed, miny, maxy, seq_id, seq_n,i_seq, org, old_regime_list)
    new_regime_list = old_regime_list;
    
    %% regime color and transparency setting
    alpha = 0.4;
    max_colors = 32;
    if max_colors > n_regime
        max_colors = n_regime;
    end
            
    %% plot segments in G regimes
    for i=0: n_regime-1
        %% regime color setting 
        % jet hsv colorcube parula
        cols = colormap(jet(max_colors));
        col=cols(mod(i, size(cols, 1)) + 1,:);

        %% load data: "segment.[i]"
        fid = fopen([datfn, num2str(i)]);
        tline = fgets(fid);
        X = [ ];
        while ischar(tline)
            dline = str2num(tline);
            if dline(3) == seq_id
                X = [X; dline];
            end
            tline = fgets(fid);
        end
        fclose(fid);

        %% if regime #i is notthing in the sequence, ignore
        if size(X,1) == 0
            continue;
        end

        %% plot segments in the regime #i
        new_regime_list(i+1) = plot_box(X, col, alpha, miny, maxy);
        plot_lines(X(:,1));
        box on
        set(gca,'YTickLabel',{''})
    end
end


function p = plot_box( X, color, alpha, miny, maxy)
    n=size(X,1);
    for i=1:n
        hold on;
        p = plot_box_aux(X(i,1), X(i,2), color, alpha, miny, maxy);
        hold off;
    end
end


function p = plot_box_aux( st, ed, color, alpha ,miny, maxy)
    btm=miny*1.1;
    top=maxy*1.1;
    p=patch([st, ed, ed, st],[btm, btm, top, top], color);
    set(p,'EdgeColor','none')
    set(p,'FaceAlpha',alpha);
    ylim([btm, top])
end


function [  ] = plot_lines( X )
    btm=-1000;
    top=1000;
    hold on;    
    for i=1: length(X)
%         plot([X(i), X(i)], [btm, top], 'black-')
    end
    %ylim([btm, top])
end


function [ Xfull ] = normalize( X )
    Xfull=X;
    for i=1: size(Xfull,2)
        X=Xfull(:,i);
        Xn=(X-min(X))/(max(X)-min(X));
        Xfull(:,i) = Xn;
    end
end


