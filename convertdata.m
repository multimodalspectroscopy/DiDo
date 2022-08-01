myDir = 'D:\uclh_data\newdata'; %gets directory
myFiles = dir(fullfile(myDir,'*.mat'));
mycount = 0;
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  if baseFileName(9:10) == 'PD'
      mycount = mycount + 1;
      fullFileName = fullfile(myDir, baseFileName);
      fprintf(1, 'Now reading %s\n', fullFileName);
      %a = dlmread(fullFileName);
      a = load(baseFileName);
      a = struct2cell(a);
      a1 = a{1,1};
      t = timetable2table(a1);
      s = size(t);
      nams = t.Properties.VariableNames;
      for m = 1:s(2)
          tab = t(:,m);
          cutname = baseFileName(1:length(baseFileName)-4);
          idv = nams{1,m}; bar = "_"; endi = ".txt";
          csvname = append(cutname, bar, idv, endi);
          writetable(tab, csvname, "Delimiter", ";");
      end
  end
  %clearvars;
  myDir = 'D:\uclh_data\newdata'; %gets directory
  myFiles = dir(fullfile(myDir,'*.mat')); 
end
mycount

