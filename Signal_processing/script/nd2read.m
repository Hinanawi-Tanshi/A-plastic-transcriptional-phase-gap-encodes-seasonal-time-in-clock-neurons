function [C1, C2] = nd2read(filename, varargin)
tic
finfo = nd2finfo(filename);
disp(['analyzing file structure used ', sprintf('%0.2f', toc), ' seconds'])

C1 = zeros(finfo.img_height, finfo.img_width, finfo.img_seq_count, 'uint16');
C2 = zeros(finfo.img_height, finfo.img_width, finfo.img_seq_count, 'uint16');

tic
Seq = finfo.file_structure(strncmp('ImageDataSeq', ...
  {finfo.file_structure(:).nameAttribute}, 12));
for frame = 1:finfo.img_seq_count
    im_ch1 = zeros(finfo.img_width, finfo.img_height, 'uint16');
    im_ch2 = zeros(finfo.img_width, finfo.img_height, 'uint16');
    fid = fopen(filename, 'r');
    fseek(fid, Seq(frame).dataStartPos, 'bof');
    
    % Image extracted from ND2 has image width defined by its first dimension.
    if finfo.padding_style == 1
      if finfo.ch_count == 2
        for ii = 1: finfo.img_height
            temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
              [finfo.ch_count finfo.img_width]);
            im_ch2(:, ii) = temp(1, :);
            im_ch1(:, ii) = temp(2, :);

            fseek(fid, 2, 'cof');
        end
      end
    else
      if finfo.ch_count == 2
        for ii = 1: finfo.img_height
            temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
              [finfo.ch_count finfo.img_width]);
            im_ch1(:, ii) = temp(1, :);
            im_ch2(:, ii) = temp(2, :);
        end 
      end
    end

    fclose(fid);

    im_ch1 = permute(im_ch1, [2 1]);
    im_ch2 = permute(im_ch2, [2 1]);
    C1(:,:,frame) = im_ch1;
    C2(:,:,frame) = im_ch2;
    
end
disp(['reading complete image data used ', sprintf('%0.2f', toc), ' seconds'])
end