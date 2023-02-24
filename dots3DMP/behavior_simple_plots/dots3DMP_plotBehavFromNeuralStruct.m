RTtask = 1;
conftask = 2;
RTCorrOnly = 1;
forTalk = 1;


data.direction = dataCell{1}.Exp.openEvents.direction';
data.coherence = dataCell{1}.Exp.openEvents.coherence';
data.scoh = data.coherence;
data.scoh(data.direction==180) = -data.scoh(data.direction==180);
data.choice = dataCell{1}.Exp.openEvents.choice'-1;
data.PDW = dataCell{1}.Exp.openEvents.pdw';
data.RT = dataCell{1}.Exp.openEvents.motionEnd'-dataCell{1}.Exp.openEvents.motionStart';
data.correct = zeros(size(data.RT));
data.correct(data.direction==0 & data.choice==1 | data.direction==180 & data.choice==0) = 1;

for N = 2:length(dataCell)
    len = length(dataCell{N}.Exp.openEvents.direction);
    data.direction(end+1:end+len) = dataCell{N}.Exp.openEvents.direction';
    data.coherence(end+1:end+len) = dataCell{N}.Exp.openEvents.coherence';
    data.choice(end+1:end+len) = dataCell{N}.Exp.openEvents.choice'-1;
    data.PDW(end+1:end+len) = dataCell{N}.Exp.openEvents.pdw';
    data.RT(end+1:end+len) = dataCell{N}.Exp.openEvents.motionEnd'-dataCell{N}.Exp.openEvents.motionStart';
end
data.scoh = data.coherence;
data.scoh(data.direction==180) = -data.scoh(data.direction==180);
data.correct = zeros(size(data.RT));
data.correct(data.direction==0 & data.choice==1 | data.direction==180 & data.choice==0) = 1;

parsedData = Dots_parseData(data,conftask,RTtask,RTCorrOnly);
cohs = unique(data.scoh);
Dots_plot(parsedData,cohs,conftask,RTtask,0,forTalk);

