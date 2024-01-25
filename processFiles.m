function processFiles(varargin)

%Validate the input
if isempty(varargin)
    selpath = uigetdir('', 'Select folder with images');

    if ~isequal(selpath, 0)
        folders = {selpath};
    else
        %Cancelled
        return;
    end

elseif ischar(varargin{1})

    folders = {varargin{1}};

elseif iscell(varargin)

    folders = varargin{1};

else

    error('Unexpected input type')

end

%Process each folder
for iDir = 1:numel(folders)

    %Find ND2 files in the folders
    files = dir(fullfile(folders{iDir}, '*.nd2'));

    %Create two CSV files to accumulate the results
    tstr = char(datetime('now'));
    tstr = strrep(tstr, ':', '-');
    tstr = strrep(tstr, ' ', '_');

    fidraw = fopen(fullfile(folders{iDir}, [tstr, '_raw.csv']), 'w');
    fprintf(fidraw, 'Filename, ObjectID, Mean GFP, Max GFP, Hit/Miss GFP, Mean TRITC, Max TRITC, Hit/Miss TRITC\n');
    
    fidsummary = fopen(fullfile(folders{iDir}, [tstr, '_summary.csv']), 'w');
    fprintf(fidsummary, 'Filename, Num Cy5 objects, Threshold GFP, Threshold TRITC, Num +GFP, Num +TRITC, Pc +GFP, Pc +TRITC\n');
   
    %Process each file
    for iFile = 1:numel(files)

        reader = BioformatsImage(fullfile(files(iFile).folder, files(iFile).name));

        try
            IGFP = getPlane(reader, 1, 'SDC-GFP', 1);
            ICy5 = getPlane(reader, 1, 'SDC-Cy5', 1);
            ITRITC = getPlane(reader, 1, 'SDC-TRITC', 1);

        catch ME

            fprintf('%s, %s. File skipped. \n', files(iFile).name, ME.message);
            continue

        end

        Mask the Cy5 channel
        maskCy5 = makeMask(ICy5);

        %Calculate raw data
        dataGFP = regionprops(maskCy5, IGFP, 'MeanIntensity', 'MaxIntensity');
        dataTRITC = regionprops(maskCy5, ITRITC, 'MeanIntensity', 'MaxIntensity');

        %Calculate thresholds
        dblIGFP = double(IGFP);
        dblITRITC = double(ITRITC);

        thGFP = mean(dblIGFP(:)) + 3 * std(dblIGFP(:));
        thTRITC = mean(dblITRITC(:)) + 3 * std(dblITRITC(:));
        % thGFP = 0.9 * max(IGFP(:));
        % thTRITC = 0.9 * max(ITRITC(:));

        %Classify each area
        posGFP = [dataGFP.MeanIntensity] > thGFP;
        posTRITC = [dataTRITC.MeanIntensity] > thTRITC;

        %Write the raw data to file
        for iObj = 1:numel(dataGFP)

            if iObj == 1
                fprintf(fidraw, '%s, %.0f, %.3f, %.3f, %.0f, %.3f, %.3f, %.0f\n', ...
                    files(iFile).name, iObj, ...
                    dataGFP(iObj).MeanIntensity, dataGFP(iObj).MaxIntensity, ...
                    posGFP(iObj), ...
                    dataTRITC(iObj).MeanIntensity, dataTRITC(iObj).MaxIntensity, ...
                    posTRITC(iObj));

            else
                fprintf(fidraw, '%s, %.0f, %.3f, %.3f, %.0f, %.3f, %.3f, %.0f\n', ...
                    ' ', iObj, ...
                    dataGFP(iObj).MeanIntensity, dataGFP(iObj).MaxIntensity, ...
                    posGFP(iObj), ...
                    dataTRITC(iObj).MeanIntensity, dataTRITC(iObj).MaxIntensity, ...
                    posTRITC(iObj));

            end

        end

        %Write summary data to file
        fprintf(fidsummary, '%s, %.0f, %.3f, %.3f, %.0f, %.0f, %.1f, %.1f\n', ...
            files(iFile).name, ...
            numel(dataGFP), thGFP, thTRITC, ...
            nnz(posGFP), nnz(posTRITC), ...
            (nnz(posGFP) / numel(dataGFP)) * 100, ...
            (nnz(posTRITC) / numel(dataTRITC)) * 100);


        %Make some output images
        % IGFPnorm = double(IGFP);
        % IGFPnorm = (IGFPnorm - min(IGFPnorm(:)))/(max(IGFPnorm(:)) - min(IGFPnorm(:)));

        ICy5norm = double(ICy5);
        ICy5norm = (ICy5norm - min(ICy5norm(:)))/(max(ICy5norm(:)) - min(ICy5norm(:)));
        Iout = imfuse(ICy5norm, bwperim(maskCy5));

        if ~exist(fullfile(files(iFile).folder, 'masked'), 'dir')
            mkdir(fullfile(files(iFile).folder, 'masked'));            
        end

        imwrite(Iout, fullfile(files(iFile).folder, 'masked', [files(iFile).name(1:end - 4), '.png']));

        % ITRITCnorm = double(ITRITC);
        % ITRITCnorm = (ITRITCnorm - min(ITRITCnorm(:)))/(max(ITRITCnorm(:)) - min(ITRITCnorm(:)));

        % Iout1 = imfuse(IGFPnorm, ICy5norm);
        % Iout1 = showoverlay(Iout1, bwperim(maskCy5), 'Color', [1 1 0]);
        % Iout1 = insertText(Iout1, [11, 12], 'GFP-Cy5', 'TextColor', 'white', 'BoxOpacity', 0);
        % 
        % Iout2 = imfuse(ITRITCnorm, ICy5norm);
        % Iout2 = showoverlay(Iout2, bwperim(maskCy5), 'Color', [1 1 0]);
        % Iout2 = insertText(Iout2, [11, 12], 'TRITC-Cy5', 'TextColor', 'white', 'BoxOpacity', 0);
        % 
        % imwrite([Iout1 Iout2], ...
        %     fullfile(files(iFile).folder, [files(iFile).name(1:end - 4), '.png']))
    end

    fclose(fidraw);
    fclose(fidsummary);

end





%
