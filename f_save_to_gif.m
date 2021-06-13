function [] = f_save_to_gif(hfig, framerate, save_gif_filename, i_frame)
% Save a figure as a frame to a gif

% Capture the plot as an image 
frame       = getframe(hfig); 
im          = frame2im(frame); 
[imind,cm]  = rgb2ind(im,256); 
% Write to the GIF File 
if i_frame == 1 
  imwrite(imind, cm, save_gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', framerate ); 
else 
  imwrite(imind, cm, save_gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', framerate ); 
end 

end

