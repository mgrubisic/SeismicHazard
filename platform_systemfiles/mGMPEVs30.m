function source = mGMPEVs30(source,Vs30)

switch source.gmpe.type
    case 'regular', str = func2str(source.gmpe.handle);
    case 'cond',    str = func2str(source.gmpe.cond.conditioning);
    case 'udm' ,    str = 'udm'; % pending
    case 'pce' ,    str = func2str(source.gmpe.handle);
end

switch str
    case 'Youngs1997' % assumptions made
        if Vs30<760 , source.gmpe.usp.media='soil';
        else        , source.gmpe.usp.media='rock';
        end
        
    case 'AtkinsonBoore2003'
        if  Vs30<180                      ,source.gmpe.usp.media='nehrpe';
        elseif and(180<=Vs30,Vs30<=360)   ,source.gmpe.usp.media='nehrpd';
        elseif and(360< Vs30,Vs30<=760)   ,source.gmpe.usp.media='nehrpc';
        elseif 760<Vs30                   ,source.gmpe.usp.media='nehrpb';
        end
        
    case 'Zhao2006'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'Mcverry2006' % assumptions made
        % used the NEHRP definition
        if     Vs30<180                  ,source.gmpe.usp.media='E';
        elseif and(180<=Vs30,Vs30< 270)  ,source.gmpe.usp.media='D';
        elseif and(270<=Vs30,Vs30< 360)  ,source.gmpe.usp.media='C';
        elseif and(360<=Vs30,Vs30<1500)  ,source.gmpe.usp.media='B';
        elseif 1500<=Vs30                ,source.gmpe.usp.media='A';
        end
        
    case 'ContrerasBoroschek2012' % assumptions made
        if Vs30<760 , source.gmpe.usp.media='soil';
        else        , source.gmpe.usp.media='rock';
        end
        
    case 'BCHydro2012'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'Arteta2018'
        if Vs30<760 , source.gmpe.usp.media='soil';
        else        , source.gmpe.usp.media='rock';
        end
        
    case 'MontalvaBastias2017'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'SiberRisk2019'
        source.gmpe.usp.Vs30=Vs30;        
        
    case 'Sadigh1997'
        if Vs30<760 , source.gmpe.usp.media='deepsoil';
        else         , source.gmpe.usp.media='rock';
        end
        
    case 'Idini2016'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'Idriss2008_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'ChiouYoungs2008_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'BooreAtkinson_2008_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'CampbellBozorgnia_2008_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'AbrahamsonSilva2008_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'AS1997h'
        if Vs30<760 , source.gmpe.usp.media='deepsoil';
        else        , source.gmpe.usp.media='rock';
        end
        
    case 'Campbell1997h'
        if     Vs30<760                  , source.gmpe.usp.media='soil';
        elseif and(760<=Vs30,Vs30<1500)  , source.gmpe.usp.media='softrock';
        end                              , source.gmpe.usp.media='hardrock';
        
    case 'I_2014_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'CY_2014_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'CB_2014_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'BSSA_2014_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'ASK_2014_nga'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'udm'
        if isfield(source.gmpe.var,'Vs30')
            source.gmpe.usp.Vs30 = Vs30;
        end
        
        if isfield(source.gmpe.var,'media')
            val = source.gmpe.var.media.value;
            val = strtrim(regexp(val,'\;','split'));
            Nlim = length(val);
            for jj=1:Nlim
                line = strtrim(regexp(val{jj},'\ ','split'));
                if eval(line{1})
                    source.gmpe.usp.media = line{2};
                end
            end
        end
        
    case 'PCE_ngawest2'
        source.gmpe.usp.Vs30=Vs30;
        
    case 'PCE_bchydro'
        source.gmpe.usp.Vs30=Vs30;
end

% These GMMs do not use site classification data

%     case 'Jaimes2006'
%     case 'Jaimes2015'
%     case 'Jaimes2016'
%     case 'GarciaJaimes2017'
