myDir = 'D:\uclh_data'; %gets directory
myFiles = dir(fullfile(myDir,'*.mat'))
mycount = 0;
bfi = [];
nir = [];
systemic = [];
variables = ['RefitDCSaDB', 'SPO2', 'HHb'];
patients = [];
sessions = [];
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  patient = baseFileName(5:7);
  session = baseFileName(25);
  
  if baseFileName(9:10) == 'PD'
      a = load(baseFileName);
      a = struct2cell(a);
      a1 = a{1,1};
      t = timetable2table(a1);
      s = size(t);
      nams = t.Properties.VariableNames;
      patients = [patients; string(patient)];
      sessions = [sessions; session];
      for item = 1:3
         
         TF = contains(nams,variables(item));
 
         if isempty(nams(TF)) == false
             if item == 1
                 bfi = [bfi ; 1];
             end
             if item == 2
                 systemic = [systemic ; 1];
             end
             if item == 3
                 nir = [nir ; 1];
             end
         else
             if item == 1
                 bfi = [bfi ; 0];
             end
             if item == 2
                 systemic = [systemic ; 0];
             end
             if item == 3
                 nir = [nir ; 0];
             end
         end
      end
  end
end

total = [patients; sessions; bfi; systemic; nir];
writematrix(total, 'D:\\uclh_data\\variables_present.csv')