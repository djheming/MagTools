classdef MagSource < matlab.mixin.Heterogeneous
    
    % MagSource: class representing a magnetic source body. This is a
    % superclass with no constructor of its own because the subclasses
    % MagBox and MagDipoles have their own particular constructors.
    %
    %
    %   Disclaimer: This code is provided as-is, has been tested only very
    %   informally, and may not always behave as intended. I find it useful for
    %   my own work, and I hope you will too, but I make no guarantees as to
    %   the accuracy or robustness of this code. This code is also actively
    %   under development and future versions may not be backward compatible.
    %
    %   Doug Hemingway (djheming@berkeley.edu)
    %   University of California, Berkeley
    %   2018-06-07
    %
    
    properties (Access=public)

        Ba = zeros(3,1) % (3x1) reference magnetic field (e.g., ambient field), for computing magnetic anomaly

    end

    methods
        
        % Display functions.
        fh = showQfieldContours( thisSource, component, varargin )
        fh = showBfieldContours( thisSource, component, varargin )
        fh = showBfieldVectors( thisSource, varargin )
        fh = showBfieldLines( thisSource, varargin )
        function drawSource( thisSource, ah, varargin )
            if isa( thisSource, 'MagBox' )
                thisSource.drawBox( ah, varargin{:} );
            elseif isa( thisSource, 'MagDipoles' )
                thisSource.drawDipoles( ah, varargin{:} );
            elseif isa( thisSource, 'MagEnsemble' )
                thisSource.drawEnsemble( ah, varargin{:} );
            end
        end

        % Setter methods.
        function thisSource = set.Ba( thisSource, newBa )
            thisSource.Ba = newBa;
            % If the user changes the ambient field for an ensemble, this
            % needs to be propagated down to the individual members of that
            % group in case any of them have induced magnetization.
            if isprop( thisSource, 'sources' )
                for k = 1 : thisSource.N
                    thisSource.sources(k).Ba = newBa;
                end
            end
        end

        % Simple data extraction functions.
        function Bcomp = getBcomponent( thisSource, B, component )
            switch component
                case 'x'
                    Bcomp = B(1,:);
                case 'y'
                    Bcomp = B(2,:);
                case 'z'
                    Bcomp = B(3,:);
                case { 't' 'total' }
                    Bcomp = sqrt( sum( B.^2, 1 ) );
                case { 'a' 'anomaly' }
                    Bcomp = sqrt(sum((B+thisSource.Ba).^2,1)) - sqrt(sum(thisSource.Ba.^2,1));
                otherwise
                    error( 'Unrecognized component type (%s)', component );
            end
        end        

    end
    
    methods (Static)

        % Basic constants.
        function result = mu_naught_over_4pi()
            result = 1e-7;
        end
        function result = mu_naught()
            result = pi * 4e-7;
        end

        % Simple translations.
        function B = B_inc_dec_to_enr( Bint, inc, dec )
            % Given a total field intensity, and a local declination
            % (rotation clockwise form north) and inclination (downward dip
            % from horizontal), compute the local vector B where the local 
            % coordinate system has +x east, +y north, +z radially upward
            % (call this the "enr" right-handed coordinate system).
            sinI = sind(inc);
            cosI = cosd(inc);
            sinD = sind(dec);
            cosD = cosd(dec);
            B = Bint * [ cosI * sinD; cosI * cosD; -sinI ];
        end

        % Display functions.
        ah = drawFieldContours( ah, survey, V, varargin )
        ah = drawBfieldVectors( ah, survey, B, varargin )
        
        % For illustration purposes.
        cmap = magcolors( n )
        
        % Test functions.
        make_figures_for_paper( figs_folder )
        
    end
    
end

