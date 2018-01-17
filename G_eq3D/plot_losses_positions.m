%to avoid plotting too prompt losses
TSTEPLOSSMIN=3000

NBTIMESTEPS=10000;
NBTIMESTAMPS=50;
ratioTS=NBTIMESTEPS/NBTIMESTAMPS;

MYPOP=find(ejected.*output.time_step_loss>TSTEPLOSSMIN);
MYPOP=MYPOP(1:1:end);
figure
hold on
grid on

for partnb=1:length(MYPOP)
    nb=MYPOP(partnb);
plot(squeeze(output.x_gc(nb,1,floor(output.time_step_loss(nb)/ratioTS)))',squeeze(output.x_gc(nb,2,floor(output.time_step_loss(nb)/ratioTS)))','.');
end