
theta_timeline_corr=theta_timeline;

end_ts=length(theta_timeline);
theta_prev=theta_timeline(1);
DTHETA_INF=-0.9*pi;
DTHETA_SUP=0.9*pi;
inf_chunk=zeros(end_ts,1);
sup_chunk=zeros(end_ts,1);

for t=2:end_ts-1
    theta_prev=theta_timeline(t-1);
    theta_next=theta_timeline(t);
    Dtheta=theta_next-theta_prev;
    if (Dtheta<DTHETA_INF)
        inf_chunk(t)=1;
    end
    if(Dtheta>DTHETA_SUP)
        sup_chunk(t)=1;
    end
end

up_indexes=find(inf_chunk);
down_indexes=find(sup_chunk);


if length(down_indexes)>0
    for t=1:length(down_indexes)
        theta_timeline_corr(down_indexes(t):end-1)=theta_timeline_corr(down_indexes(t):end-1)-2*pi;
    end
end

if length(up_indexes)>0
    for t=1:length(up_indexes)
        theta_timeline_corr(up_indexes(t):end-1)=theta_timeline_corr(up_indexes(t):end-1)+2*pi;
    end
end

Dtheta=theta_timeline(end)-theta_timeline(end-1);
theta_timeline_corr(end)=theta_timeline_corr(end-1)+Dtheta;
