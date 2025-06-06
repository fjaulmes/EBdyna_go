function Beam=read_TRANSP_namelist(file)

%read_namelist should be in the path
%test: file='/compass/home/zadvitskiy/matlab/NUR/99999I38TR.DAT';
%test: file='/compass/home/zadvitskiy/matlab/NUR/test.DAT';
%test: addpath('/compass/home/zadvitskiy/matlab/rw_namelist/');

Beam=struct;

indata=read_namelist_TRANSP_NUR(file);

try
    Beam.nbeam=indata.TRANSP.NBEAM;
catch
    indata.TRANSP.NBEAM=1;
    Beam.nbeam=indata.TRANSP.NBEAM;
end

for i=1:indata.TRANSP.NBEAM
    %grid
    try
        if numel(indata.TRANSP.BMWIDR)>1
            Beam.bmwid.r(i)=indata.TRANSP.BMWIDR(i)*1e-2;  %source/(last grid) half width [m]
            Beam.bmwid.z(i)=indata.TRANSP.BMWIDZ(i)*1e-2;  %source/(last grid) half width [m]
        else
            Beam.bmwid.r(i)=indata.TRANSP.BMWIDR*1e-2;  %source/(last grid) half width [m]
            Beam.bmwid.z(i)=indata.TRANSP.BMWIDZ*1e-2;  %source/(last grid) half width [m]
        end
    catch
        warning('source/(last grid) half width is not defined in the namelist');
    end
    %shape option
    try
        if numel(indata.TRANSP.NBSHAP)>1
            Beam.bmwid.opt(i)=indata.TRANSP.NBSHAP(i);  %option shape 0-circ, 1-rect
        else
            Beam.bmwid.opt(i)=indata.TRANSP.NBSHAP;  %option shape 0-circ, 1-rect
        end
    catch
        warning('source/(last grid) half width is not defined in the namelist');
    end
    %focal distance
    try
        if numel(indata.TRANSP.FOCLR)>1
            Beam.focl.r(i)=indata.TRANSP.FOCLR(i)*1e-2;  %distance from the source to optical focual point [m]
            Beam.focl.z(i)=indata.TRANSP.FOCLZ(i)*1e-2;  %distance from the source to optical focual point [m]
        else
            Beam.focl.r(i)=indata.TRANSP.FOCLR*1e-2;  %distance from the source to optical focual point [m]
            Beam.focl.z(i)=indata.TRANSP.FOCLZ*1e-2;  %distance from the source to optical focual point [m]
        end
    catch
        warning('focal distance is not defined in the namelist');
    end
    %source divergence
    try
        if numel(indata.TRANSP.DIVR)>1
            Beam.div.r(i)=indata.TRANSP.DIVR(i); %source divergence [rad]
            Beam.div.z(i)=indata.TRANSP.DIVZ(i); %source divergence [rad]
        else
            Beam.div.r(i)=indata.TRANSP.DIVR; %source divergence [rad]
            Beam.div.z(i)=indata.TRANSP.DIVZ; %source divergence [rad]
        end
    catch
        warning('source divergence is not defined in the namelist');
    end
    %distance from the source to the apperture
    try
        if numel(indata.TRANSP.XLBAPA)>1
            Beam.app.l(i)=indata.TRANSP.XLBAPA(i)*1e-2; %distance from the source to the apperture [m]
        else
            Beam.app.l(i)=indata.TRANSP.XLBAPA*1e-2; %distance from the source to the apperture [m]
        end
    catch
        warning('distance from the source to the apperture is not defined in the namelist');
    end
    %apperture half width
    try
        if numel(indata.TRANSP.BMWIDR)>1
            Beam.app.r(i)=indata.TRANSP.BMWIDR(i)*1e-2; %apperture half width [m]
            Beam.app.z(i)=indata.TRANSP.BMWIDZ(i)*1e-2; %apperture half width [m]
        else
            Beam.app.r(i)=indata.TRANSP.BMWIDR*1e-2; %apperture half width [m]
            Beam.app.z(i)=indata.TRANSP.BMWIDZ*1e-2; %apperture half width [m]
        end
    catch
        warning('apperture half width is not defined in the namelist');
    end
    %Tangency radius [m]
    try
        if numel(indata.TRANSP.RTCENA)>1
            Beam.R_t(i)=indata.TRANSP.RTCENA(i)*1e-2; %Tangency radius [m]
        else
            Beam.R_t(i)=indata.TRANSP.RTCENA*1e-2; %Tangency radius [m]
        end
    catch
        warning('Tangency radius is not defined in the namelist');
    end
    %distance from the source to the tangency point [m]
    try
        if numel(indata.TRANSP.XLBAPA)>1
            Beam.d_s(i)=indata.TRANSP.XLBTNA(i)*1e-2;% distance from the source to the tangency point [m]
        else
            Beam.d_s(i)=indata.TRANSP.XLBTNA*1e-2;% distance from the source to the tangency point [m]
        end
    catch
        warning('Tangency radius is not defined in the namelist');
    end
    %source elevation over midplane [m] and 1st aperture elevation over midplane [m]
    try
        if numel(indata.TRANSP.XYBSCA)>1
            Beam.Z_s(i)=indata.TRANSP.XYBSCA(i)*1e-2; %source elevation over midplane [m]
            Beam.Z_a(i)=indata.TRANSP.XYBAPA(i)*1e-2; %1st aperture elevation over midplane [m]
        else
            Beam.Z_s(i)=indata.TRANSP.XYBSCA*1e-2; %source elevation over midplane [m]
            Beam.Z_a(i)=indata.TRANSP.XYBAPA*1e-2; %1st aperture elevation over midplane [m]
        end
    catch
        warning('source elevation over midplane and 1st aperture elevation over midplane are not defined in the namelist');
    end
    %beam speace atomic number
    try
        if numel(indata.TRANSP.ABEAMA)>1
            Beam.Ab(i)=indata.TRANSP.ABEAMA(i); %beam speace atomic number
        else
            Beam.Ab(i)=indata.TRANSP.ABEAMA; %beam speace atomic number
        end
    catch
        warning('beam speace atomic number is not defined in the namelist');
    end
    %Energy
    try
        if numel(indata.TRANSP.EINJA)>1
            Beam.E(i)=indata.TRANSP.EINJA(i); %Energy
        else
            Beam.E(i)=indata.TRANSP.EINJA; %Energy
        end
    catch
        warning('beam energy is not defined in the namelist');
    end
    %1/2 E 1/3 E current fraction
    try
        if numel(indata.TRANSP.FHALFA)>1
            Beam.frac_E1(i)=indata.TRANSP.FFULLA(i); %current fraction
            Beam.frac_E2(i)=indata.TRANSP.FHALFA(i); %current fraction
        else
            Beam.frac_E1(i)=indata.TRANSP.FFULLA; %current fraction
            Beam.frac_E2(i)=indata.TRANSP.FHALFA; %current fraction
        end
    catch
        warning('1/2 E 1/3 E current fraction are not defined in the namelist');
    end
    %power
    try
        if numel(indata.TRANSP.PINJA)>1
            Beam.P(i)=indata.TRANSP.PINJA(i); %current fraction
            
        else
            Beam.P(i)=indata.TRANSP.PINJA; %current fraction
            
        end
    catch
        warning('power not defined in the namelist');
    end
    %toroidal angle
    try
        if numel(indata.TRANSP.XBZETA)>1
            Beam.phi(i)=indata.TRANSP.XBZETA(i); %toroidal angle            
        else
            Beam.phi(i)=indata.TRANSP.XBZETA; %toroidal angle           
        end
    catch
        warning('toroidal angle is defined in the namelist');
    end
end

return
end


