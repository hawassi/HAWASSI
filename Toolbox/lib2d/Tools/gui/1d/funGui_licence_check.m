function IdLicence=licence_check(handles,pathLic)
load(pathLic);
Regkey=deserialize(regkey);
startdate=Regkey.startdate;

RegPassword=Regkey.password{1,1};

startdate = datenum(startdate);
todaysdate = datenum(date);

global IDdemo

if IDdemo==0
    Licence_List{1}='1U56-1BC3-0UVD-KJM1';
    Licence_List{2}='1STV-1LMF-056D-54B2';
    Licence_List{3}='16BU-1AB9-08CA-A0C3';
    Licence_List{4}='1BAA-154A-052Z-5BZ1';
    Licence_List{5}='1Z69-1Z3T-0KM1-KUT2';
    Licence_List{6}='1CBK-1RUK-0ER9-7BA3';
else
    Licence_List{1}='1BJK-087A-0UQ9-KKU1';
    Licence_List{2}='156Q-0LU7-0BK4-KUA1';
    Licence_List{3}='1KU9-0HFI-07K9-2WT1';
    Licence_List{4}='1LMI-0NEX-0978-KKL2';
    Licence_List{5}='1KYU-0J67-0KQA-LTU3';
end
           

[~,result] = dos('getmac');
mac = result(160:176);


if strcmp(mac,Regkey.mac)==1 
    if Regkey.flag==0
        if strcmp(RegPassword,Licence_List{1})==1
            Regkey.duration=365;
            Regkey.flag    =1;
            set(handles.Activation,'userdata',1);
            regkey=serialize(Regkey);
            save(pathLic,'regkey')
            IdLicence=1;
        elseif strcmp(RegPassword,Licence_List{2})==1
            Regkey.duration=365;
            Regkey.flag    =1;
            set(handles.Activation,'userdata',1);
            regkey=serialize(Regkey);
            save(pathLic,'regkey')
            IdLicence=1;
        elseif strcmp(RegPassword,Licence_List{3})==1
            Regkey.duration=365;
            Regkey.flag    =1;
            set(handles.Activation,'userdata',1);
            regkey=serialize(Regkey);
            save(pathLic,'regkey')
            IdLicence=1;
        elseif strcmp(RegPassword,Licence_List{4})==1
            Regkey.duration=30;
            Regkey.flag    =1;
            set(handles.Activation,'userdata',1);
            regkey=serialize(Regkey);
            save(pathLic,'regkey')
            IdLicence=1;
        elseif strcmp(RegPassword,Licence_List{5})==1
            Regkey.duration=30;
            Regkey.flag    =1;
            set(handles.Activation,'userdata',1);
            regkey=serialize(Regkey);
            save(pathLic,'regkey')
            IdLicence=1;
        elseif strcmp(RegPassword,Licence_List{6})==1
            Regkey.duration=30;
            Regkey.flag    =1;
            set(handles.Activation,'userdata',1);
            regkey=serialize(Regkey);
            save(pathLic,'regkey')
            IdLicence=1;   
        else
            h = errordlg('Licence is invalid!');
            set(h,'units','normalized','Position',[0.7  0.4 0.15 0.1]);
            IdLicence=0;
            set(handles.Activation,'userdata',0);
%             delete(pathLic);
        end
        
    elseif Regkey.flag==1;
       if any(strcmp(RegPassword,Licence_List)==1) 
        DurationLicence=Regkey.duration;
        if (todaysdate - startdate >= DurationLicence)
            IdLicence=0;
             h = errordlg('Licence has expired!');
               set(h,'units','normalized','Position',[0.7  0.4 0.15 0.1]);
        else
          if todaysdate-startdate>DurationLicence-10
                daysleft=DurationLicence-(todaysdate-startdate);
                if daysleft>1
                h = warndlg(['Licence expires in ',num2str(daysleft),' days']);
                else
                h = warndlg(['Licence expires in ',num2str(daysleft),' day']);    
                end
                set(h,'units','normalized','Position',[0.7  0.4 0.15 0.1]);
          end
            
            IdLicence=1;  
        end
       else
           
            h = errordlg('Licence is invalid!'); %%for update a new software
            set(h,'units','normalized','Position',[0.7  0.4 0.15 0.1]);
            IdLicence=0;
            set(handles.Activation,'userdata',0);   
       end
    end
else

h = errordlg('Licence is invalid!');
set(h,'units','normalized','Position',[0.7  0.4 0.15 0.1]);
    IdLicence=0;
end
