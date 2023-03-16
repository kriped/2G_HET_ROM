function line = get_axial(mode)
[~,max_idx] = max(mode,[],'all','linear');
    [x_max,y_max,~]=ind2sub(size(mode),max_idx);
    line(:) = mode(x_max,y_max,:);
end