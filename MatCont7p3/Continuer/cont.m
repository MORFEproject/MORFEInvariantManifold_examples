function [xout, vout, sout, hout, fout] = cont(varargin)
%
% CONTINUE(cds.curve, x0, v0, options)
%
% Continues the curve from x0 with optional directional vector v0
% options is a option-vector created with CONTSET
% The first two parameters are mandatory.

% WM: Rearranged many things throughout the code to put it
% all in a more logical order
global cds sOutput
%    Interactive plotting
gui_enabled = ~isempty(sOutput) && sOutput.enabled;

if gui_enabled
    plotpoints = sOutput.getPlotPointsInterval();
    sOutput.setStatus('Computing ...');
    drawnow;
else
    plotpoints = 0;
    sOutput = struct();
    sOutput.enabled = false;
    sOutput.setStatus = @(varargin) 0;  %dummy function, does nothing, takes any variables
    sOutput.endRun = @(varargin) 0;
    sOutput.output = @(varargin) 0;
    sOutput.checkPauseResumeStop = @(varargin) 0;
end


if (nargin < 6)% case not extend
    if isfield(cds,'symjac')
        tmpjacs = [cds.symjac cds.symhess];
    end
    cds = [];
    cds.curve = '';

    lastwarn('');
    
    % Parse given options
    [cds.curve, x0, v0, opt] = ParseCommandLine(varargin{:});
    
    if size(x0,2) ~= 1
        x0 = x0';
    end
    cds.ndim = length(x0);
    curvehandles = feval(cds.curve);
    cds.curve_func = curvehandles{1};
    cds.curve_defaultprocessor = curvehandles{2};
    cds.curve_options = curvehandles{3};
    cds.curve_jacobian = curvehandles{4};
    cds.curve_hessians = curvehandles{5};
    cds.curve_testf = curvehandles{6};
    cds.curve_userf = curvehandles{7};
    cds.curve_process = curvehandles{8};
    cds.curve_singmat = curvehandles{9};
    cds.curve_locate = curvehandles{10};
    cds.curve_init = curvehandles{11};
    cds.curve_done = curvehandles{12};
    cds.curve_adapt = curvehandles{13};
    % Read out all options or set default values
    %
    % Get options from curve file
    options = feval(cds.curve_options);
%     options = feval(cds.curve, 'options');
    % Merge options from curve with cmdline, cmdline overrides curve
    cds.options = contmrg(opt, options);
    
    cds.options.MoorePenrose      = contget(cds.options, 'MoorePenrose', 1);
    cds.options.SymDerivative     = contget(cds.options, 'SymDerivative', 0);
    cds.options.SymDerivativeP    = contget(cds.options, 'SymDerivativeP', 0);
    cds.options.Increment         = contget(cds.options, 'Increment', 1e-5);
    cds.options.MaxNewtonIters    = contget(cds.options, 'MaxNewtonIters', 3);
    cds.options.MaxCorrIters      = contget(cds.options, 'MaxCorrIters', 10);
    cds.options.MaxTestIters      = contget(cds.options, 'MaxTestIters', 10);
    cds.options.FunTolerance      = contget(cds.options, 'FunTolerance', 1e-6);
    cds.options.VarTolerance      = contget(cds.options, 'VarTolerance', 1e-6);
    cds.options.TestTolerance     = contget(cds.options, 'TestTolerance', 1e-5);
    cds.options.TSearchOrder      = contget(cds.options, 'TSearchOrder', 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Userfunctions         = contget(cds.options, 'Userfunctions',0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Singularities         = contget(cds.options, 'Singularities', 0);
    WorkSpace             = contget(cds.options, 'WorkSpace', 0);
    Backward              = contget(cds.options, 'Backward',0);
    CheckClosed           = contget(cds.options, 'CheckClosed', 50);
    npoints               = contget(cds.options, 'MaxNumPoints', 300);
    Adapt                 = contget(cds.options, 'Adapt',3);
    
    Locators              = contget(cds.options, 'Locators', []);
    IgnoreSings           = contget(cds.options, 'IgnoreSingularity', []);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    UserInfo             = contget(cds.options, 'UserfunctionsInfo', []);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cds.h                 = contget(cds.options, 'InitStepsize' , 0.01);
    cds.h_max             = contget(cds.options, 'MaxStepsize'  , 0.1);
    cds.h_min             = contget(cds.options, 'MinStepsize'  , 1e-5); 
    
    cds.oldJacFull = [];
    cds.oldJac = [];
    cds.oldJacX = [];
    
    if exist('tmpjacs','var')
        cds.symjac = tmpjacs(1);
        cds.symhess = tmpjacs(2);
    else
        cds.symjac = 0;
        cds.symhess = 0;
    end
    
    % Set algorithm specific variables & declarations
    %
    cds.h_inc_fac = 1.3;
    cds.h_dec_fac = 0.5;
    cds.dir_check_step = 3;
    cds.dir_check_angle = 0.9;
    
    x1   = zeros(cds.ndim,1);           %
    x2   = zeros(cds.ndim,1);           % convenient to have these 4 vectors,
    v1   = zeros(cds.ndim,1);           % x1 and x2 are last two points
    v2   = zeros(cds.ndim,1);           % v1 and v2 the corresponding tangent vectors
    
    i    = 0;                       % our iteration parameter
    ind  = [];                       % our iteration parameter to plot
    it   = 0;                       % number of newton steps
    
    % Get info about testfunctions and singularities
    %
    
    cds.nActTest = 0;
    if Singularities
        [cds.S,SingLables] = feval(cds.curve_singmat);
        [cds.nSing, cds.nTest] = size(cds.S);
        Singularities = cds.nSing>0 & cds.nTest>0;
        % setup testfunction variables and stuff
        %
        
        if Singularities
            % Ignore Singularities
            cds.S(IgnoreSings,:) = 8;
            ActSing = setdiff(1:cds.nSing, IgnoreSings);
            cds.ActSing = ActSing;
            cds.nActSing = length(ActSing);
            
            % Which test functions must vanish somewhere in S?
            ActTest = find( sum((cds.S==0),1) > 0 );
            cds.ActTest = ActTest;
            cds.nActTest = length(ActTest);
            % WM: Build matrix with indices of testfunctions which
            % should be zero for each singularity            
            cds.SZ = zeros(cds.nTest+1,cds.nActSing+1);
            ml = 2;
            for si=ActSing(1:cds.nActSing)
                t = find( cds.S(si,:) == 0 )';
                l = size(t,1);
                cds.SZ(1:l,si) = t;
                ml = max(ml,l);
            end
            cds.SZ = cds.SZ(1:ml,:);
             cds.atv = 1;    cds.testvals = zeros(2, cds.nActTest);  
            % 1st row: testvals at x1, 2nd: testvals at x2
            cds.testzero = zeros(cds.ndim, cds.nActTest);  % coords where testf is zero
            cds.testvzero = zeros(cds.ndim, cds.nActTest);  % v where testf is zero
            
            if isempty(ActTest)
                Singularities = 0;
            end
        end
    end
    if Userfunctions
        cds.nUserf = size(UserInfo,2);cds.utv = 1;
        cds.uservals = zeros(2,cds.nUserf);%1st row testvals at x1,2nd testvals at x2
    end
    
    
    xout = zeros(cds.ndim,npoints);     % each point on curve
    vout = zeros(cds.ndim,npoints);     % all tangent vectors
    %hout = zeros(2+cds.nActTest+size(UserInfo,2),npoints); % keep track of all stepsizes & testfunctions
    sout = [];                          % all occured singularities
%     fout = [];
    
    %interactive plotting
 
    %if (sOutput.enabled)
    %    sOutput.output(xout,[],[],[],1);
    %end

   
    % Algorithm starts here
    %

    StartTime = clock;
        
    % restarting point...
    
    % Init user space
    %
    if WorkSpace
        if feval(cds.curve_init, x0, v0)~=0
            sOutput.setStatus('Error !!!');
            errordlg('Initializer failed.');
            sOutput.endRun();
            return;
        end
    end
    % if x0 and v0 known: assume point on curve
    % if v0 unknown: correct x0

    if isempty(v0)
        [x0, v0] = CorrectStartPoint(x0, v0);
        if isempty(x0)
            printconsole('No convergence at x0.\n') 
            printconsole('elapsed time  = %.1f secs\n', etime(clock, StartTime));
            printconsole('0 npoints\n');            

            errordlg('No convergence at x0.');   
            xout = []; vout = []; sout = []; hout = []; fout = [];
            sOutput.setStatus('No convergence at x0');
            sOutput.endRun();
            return;            
        end
    end
    if Backward
        v0 = -v0;
    end
    
    % WM: added call to the default processor
    s.index = 1;
    s.label = '00';
    s.data  = [];
    s.msg   = 'This is the first point of the curve';
    
    if gui_enabled; s.internallabel = []; end
    
    [failed,f,s] = DefaultProcessor(x0,v0,s);
    printconsole('first point found\n');%printv(x0);
    printconsole('tangent vector to first point found\n'); %printv(v0);
    x1 = x0;
    v1 = v0;
    cds.testvals = [];
    cds.uservals = [];
    if Singularities
        % WM: calculate all testfunctions at once
        [tfvals,failed] = EvalTestFunc(ActTest,x0,v0);
        cds.testvals(2,:) = tfvals(ActTest);
        if ~isempty(failed)
            sOutput.setStatus('Error');
            printconsole('Evaluation of test functions failed at starting point.\n');
            errordlg('Evaluation of test functions failed at starting point.');
            sOutput.endRun();
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Userfunctions
        [ufvals,failed]=feval(cds.curve_userf, UserInfo,1:cds.nUserf, x0, v0);
        cds.uservals(2,:)=ufvals;
        if ~isempty(failed)
            sOutput.setStatus('Error');
            printconsole('Evaluation of userfunctions failed at starting point.\n');
            errordlg('Evaluation of userfunctions failed at starting point.');
            sOutput.endRun();
        end        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xout(:,1) = x0;
    vout(:,1) = v0;
    sout = [s]; 
    if Singularities & ~Userfunctions
        hout(:,1) = [0;0;cds.testvals(2,:)'];
    elseif Userfunctions & ~Singularities
        hout(:,1) = [0;0;cds.uservals(2,:)'];
    elseif Userfunctions & Singularities
        hout(:,1) = [0;0;cds.uservals(2,:)';cds.testvals(2,:)'];
    else
        hout(:,1) = [0;0];
    end
    sf = size(f,1);
    if sf > 0
        fout = zeros(size(f,1),npoints);
    else
        fout = zeros(1,npoints);
    end
    fout(:,1) = f;
    i   = 1;
    ind = 1;
else
    xout = varargin{1};  i  = size(xout,2);
    x0 = xout(:,1);   x1 = xout(:,end);
    vout = varargin{2};
    v1 = vout(:,end);    ind = 1;
    npoints = contget(cds.options, 'MaxNumPoints', 300);
    xout = [ xout zeros(cds.ndim,npoints)];
    vout = [ vout zeros(cds.ndim,npoints)]; 
    sout =  varargin{3};
    if gui_enabled
        sout(1).internallabel = [];  %add field if not already in struct.
    end
    hout =  varargin{4};
    fout = [ varargin{5} zeros(size(varargin{5},1),npoints)];
    cds  = varargin{6};
    npoints = npoints + i;
    sout(end) = [];
    s = sout(end);
    
  
    Userfunctions         = contget(cds.options, 'Userfunctions',0);
    Singularities         = contget(cds.options, 'Singularities', 0);
    WorkSpace             = contget(cds.options, 'WorkSpace', 0);
    CheckClosed           = contget(cds.options, 'CheckClosed', 50);
    Adapt                 = contget(cds.options, 'Adapt', 1);
    
    Locators              = contget(cds.options, 'Locators', []);
    IgnoreSings           = contget(cds.options, 'IgnoreSingularity', []);
    UserInfo              = contget(cds.options, 'UserfunctionsInfo', []);
  
    if Singularities
        ActTest = cds.ActTest;
        ActSing = cds.ActSing;
        [cds.S,SingLables] = feval(cds.curve_singmat);
        
        if isempty(ActTest)
            Singularities = 0;
        end
    end
    if Userfunctions
        cds.nUserf = size(UserInfo,2);cds.utv = 1;
        cds.uservals = zeros(2,cds.nUserf);%1st row testvals at x1,2nd testvals at x2
    end
    if Userfunctions
        [ufvals,failed]=feval(cds.curve_userf, UserInfo,1:cds.nUserf, x1, v1);
        cds.uservals(2,:)=ufvals;
        if ~isempty(failed)
            sOutput.setStatus('Error');
            printconsole('Evaluation of userfunctions failed at starting point.\n'); 
            errordlg('Evaluation of userfunctions failed at starting point.'); 
            sOutput.endRun();

        end        
    end


    
    %{
    if (sOutput.enabled)
        sOutput.output(xout,[],[],[],1);
    end
    %}
    
    StartTime = clock;
    printconsole('start computing extended curve\n');
end
while i < npoints
    % correction
    %
    corrections = 1;
    while 1
      % predict next point  
      xpre = x1 + cds.h * v1;
      
      [x2, v2, it] = newtcorr(xpre, v1);
        
      if ~isempty(x2) && ((i < cds.dir_check_step) || (v1'*v2 > cds.dir_check_angle))
            % we may have found a new point, call default processor
            % to initialize some things for test functions
            [failed,f] = DefaultProcessor(x2,v2);
            if ~failed
                if Singularities
                    % WM: evaluate all testfunctions at the same time
                    [tfvals,failed] = EvalTestFunc(ActTest,x2,v2);
                    cds.testvals(cds.atv,:) = tfvals(ActTest);
                end
                % WM: if all done succesfully then we have our new point
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                if Userfunctions     
                    % evaluate all userfunctions at the same time
                    [ufvals,failed]=feval(cds.curve_userf, UserInfo,1:cds.nUserf, x2, v2);
                    cds.uservals(cds.utv,:)=ufvals;
                end
                if isempty(failed)|~failed
                    break
                end
            end  
      end
      % decrease stepsize if allowed
      if cds.h > cds.h_min
          cds.h = max(cds.h_min, cds.h*cds.h_dec_fac);
          corrections = corrections + 1;
      else      % if not then fail
          
          printconsole('Current step size too small (point %d)\n',i);
          xout = xout(:,1:i);
          vout = vout(:,1:i);
          hout = hout(:,1:i);
          fout = fout(:,1:i);
          % Last point is also singular
          % WM: Added call to default processor
          s.index = i;
          s.label = '99';
          %s.msg   = 'This is the last point on the curve';
          s.data  = [];
          s.msg   = 'This is the last point on the curve';
          if gui_enabled; s.internallabel = []; end
          [failed,f,s] = DefaultProcessor(xout(:,i), vout(:,i), s);

          
          sout = [sout; s];
          fout = [fout  f];
          EndTime = clock;
          printconsole('elapsed time  = %.1f secs\n', etime(EndTime, StartTime));
          printconsole('npoints curve = %d\n', i);
          string = sprintf('%.1f secs\n', etime(EndTime, StartTime));
          
          %sOutput  'string' duration
          sOutput.setStatus('Current step size too small');
          sOutput.endRun();
          return;
      end
   end
  % Singularities
  %
  if Singularities
      % WM: the testvals arrays are not copied anymore, instead
      % after every iteration the function of both is swapped
      cds.atv = 3-cds.atv;
      % WM: use sign function and compare instead of multiply (for speed).
      testchanges = (sign(cds.testvals(1,:)) ~= sign(cds.testvals(2,:)));
      if any(testchanges)
          % Change by WG     
           testidx = [ActTest(find( testchanges )) 0]';
           % check if sing occured
           % WM: detect all singularities with a single ismember() call
%         stz = ismember(cds.SZ,testidx);
%           singsdetected(ActSing) = all(stz(:,ActSing));
      % DSB: Singularity is detected if
      %  - Every crossing that is required occurs
      %  - Every crossing that is not required does not occur
      %
       S_true  = +(cds.S(:,ActTest)' == 0);  % Required crossings matrix
       S_false = +(cds.S(:,ActTest)' == 1);  % Required noncrossings matrix
       all_sings_detected = (testchanges*S_true ==sum(S_true))&(~testchanges*S_false == sum(S_false));
       singsdetected(ActSing) = all_sings_detected(ActSing);
          if any(singsdetected)
              % singularity detected!
              
              singsdetected = find(singsdetected==1);
              cds.testzero = zeros(cds.ndim,cds.nTest);
              cds.testvzero = zeros(cds.ndim,cds.nTest);
              % locate zeros of all testf which changed sign
              %
              xss = [];  % x of singularites
              vss = [] ; % v of idem
              testfound = [];% indices of found zeros of test functions
              
              sid = [];  % id of idem
              for si=singsdetected
                  %         printconsole('Singularity %d detected ... ', si);
                  if ismember(si, find(Locators==1))   % do we have a locator?
%                       printconsole('using user locator\n');
                      [xs,vs] = feval(cds.curve_locate, si, x1, v1, x2, v2);
                      lit=0;
                      if ~isempty(xs) &(norm(xs-(x1+x2)/2)<2*norm(x2-x1))
                          xss = [xss xs];
                          vss = [vss vs];
                          sid = [sid si];
                      end
                  else             
                      % locate zeros of test functions if not already computed
                      for ti=1:cds.nTest   
                          if cds.S(si,ti)==0 & ~ismember(ti, testfound) 
                              [xtf,vtf,lit] = LocateTestFunction(ti, x1, v1, x2, v2);
                              
                              if ~isempty(xtf)&(norm(xtf-(x1+x2)/2)<2*norm(x2-x1))                               
                                  cds.testzero(:,ti) = xtf;
                                  cds.testvzero(:,ti) = vtf;
                                  testfound = [testfound ti];
                      
%                               else
%                                   printconsole('A testfunction for Singularity %d failed to converge ... \n', si);
                              end
                          end
                      end
                      % now we have all zeros/testfunctions we need for
                      % this singularity
                      if any(ismember(testfound, find(cds.S(si,:)==0)))
                          [xs,vs] = LocateSingularity(si);
                      else
                          xs = [];
                          vs = [];
                      end
                      if ~isempty(xs)
                          xss = [xss xs];
                          vss = [vss vs];
                          sid = [sid si];
                      end                          
                  end
              end %end of detect/locate loop   
              if ~isempty(sid)         % sort
                  [xss,vss,sid] = xssort(x1, xss, vss, sid);
                  % WM: moved out of loop for speed
                  sids = 1:length(sid);
                  isids = i+sids;
                  xout(:,isids) = xss(:,sids);
                  vout(:,isids) = vss(:,sids);
                  % WM: moved call to default processor up before
                  % the special processor
                  for si=sids                      
                      i = i+1; ind = [ind i];s = [];
                      s.index = i;
                      s.label = SingLables(sid(si),:);
                      if gui_enabled; s.internallabel = []; end
                      [failed,sf,s] = DefaultProcessor(xss(:,si), vss(:,si), s);
                      [tfvals,failed] = EvalTestFunc(ActTest,xss(:,si),vss(:,si));
                      s.data.testfunctions = tfvals(ActTest);
                      [failed,s] = feval(cds.curve_process, sid(si), xss(:,si), vss(:,si), s );
                      if Userfunctions
                          [ufvals,failed] =feval(cds.curve_userf, UserInfo, 1:cds.nUserf, xss(:,si), vss(:,si));
                          s.data.userfunctions=ufvals;
%                           hout(:,i)=[cds.h;0;s.data.userfunctions';s.data.testfunctions'];
                          hout(:,i)=[0;lit;s.data.userfunctions';s.data.testfunctions'];
                      else
                          % XXX lit = # localisationsteps XXX
%                           hout(:,i) = [cds.h;lit;s.data.testfunctions'];
                          hout(:,i) = [0;lit;s.data.testfunctions'];
                      end
                      sout = [sout; s];
                      fout(:,i) = sf;
                      sOutput.output(xout,s,hout,fout,ind);
                      sOutput.setStatus(s.msg);
                      drawnow;
                      doStop = sOutput.checkPauseResumeStop(s.index , s.msg ,i);
                      if doStop~=0, npoints = i;end  %npoints = i ---> end of mainloop
                      
                  end
              end % end of loop over singularities
          end
      end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if Userfunctions
      cds.utv = 3-cds.utv;
      % WM: use sign function and compare instead of multiply (for speed).
      userchanges = sign(cds.uservals(1,:)) ~= sign(cds.uservals(2,:));
      if any(userchanges)
          % Change by WG     
          useridx = find(userchanges);
          % cds.SZ=[2;0]
          % check if sing occured
          % WM: detect all singularities with a single ismember() call
          if any(useridx)
              % singularity detected!
              %   cds.userzero = zeros(cds.ndim,size(useridx,2));
              %   cds.uservzero = zeros(cds.ndim,size(useridx,2));
              % locate zeros of all testf which changed sign%
              xus = [];  % x of singularites
              vus = []; % v of idem
              uid = [];  % id of idem
              % locate zeros of userfunctions if not already computed
              for ti=1:size(useridx,2)
                  [xtf,vtf,lit] = LocateUserFunction(UserInfo(useridx(ti)),useridx(ti), x1, v1, x2, v2);
                  if ~isempty(xtf)
                      xus = [xus xtf];
                      vus = [vus vtf];
                      uid = [uid useridx(ti)];          
                  end  
              end
          end %end of detect/locate loop
          if ~isempty(uid)         % sort
              [xus,vus,uid] = xssort(x1, xus, vus, uid);
              % WM: moved out of loop for speed
              uids = 1:length(uid);
              uuids = i+uids;
              xout(:,uuids) = xus(:,uids);
              vout(:,uuids) = vus(:,uids);
              % WM: moved call to default processor up before
              % the special processor
              for ui=uids
                  i = i+1; ind = [ind,i]; s = [];                         
                  s.index = i;
                  s.label = UserInfo(uid(ui)).label;
                  if gui_enabled; s.internallabel = []; end
                  [failed,uf,s] = DefaultProcessor(xus(:,ui), vus(:,ui), s);
                  [ufvals,failed] =feval(cds.curve_userf, UserInfo,1:cds.nUserf, xus(:,ui), vus(:,ui));
                  s.data.userfunctions = ufvals;
                  if Singularities
                      [tfvals,failed] = EvalTestFunc(ActTest,xus(:,ui),vus(:,ui));
                      s.data.testfunctions = tfvals(ActTest);
                      hout(:,i) = [0;lit;s.data.userfunctions';s.data.testfunctions'];
                  else
                      hout(:,i) = [0;lit;s.data.userfunctions'];
                  end
                  printconsole('label = %s, x = ', s.label); printv((xus(:,ui)));
                  s.msg  = sprintf('%s',UserInfo(uid(ui)).name);
                  sout = [sout; s];
                  fout(:,i) = uf;
                  if (sOutput.enabled)
                      sOutput.output(xout,s,hout,fout,ind);
                      doStop = sOutput.checkPauseResumeStop(s.index , s.msg ,i);
                      if doStop~=0, npoints = i;end  %npoints = i ---> end of mainloop
                      ind =[];
                  end
                  sOutput.setStatus(s.msg);
              end
          end
      end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (mod(i,Adapt)==0)
   % WM: Adapt can now return new x and v as well
    [res,x2,v2] = feval(cds.curve_adapt, x2, v2);
    [failed,f] = DefaultProcessor(x2,v2);
    if res==1 & Singularities
      % recompute testvals
      [tfvals,failed] = EvalTestFunc(ActTest,x2,v2);
  	  cds.testvals(3-cds.atv,:) = tfvals(ActTest);
    end
  end
  
  % WM: Moved this up a bit the prevent several calculations of i+1
  i   = i+1;
  ind = [ind,i];
  xout(:,i) = x2;
  vout(:,i) = v2;
  %%%%%%%%%%%%%%%%%%%
  if Singularities & ~Userfunctions
      hout(:,i) = [cds.h;it;cds.testvals(3-cds.atv,:)'];
  elseif Userfunctions & ~Singularities
      hout(:,i) = [cds.h;it;cds.uservals(3-cds.utv,:)'];
  elseif Userfunctions & Singularities
      hout(:,i) = [cds.h;it;cds.uservals(3-cds.utv,:)';cds.testvals(3-cds.atv,:)'];
  else
      hout(:,i) = [cds.h;it];
  end
  %%%%%%%%%%%%%%%%%%%
  fout(:,i) = f;
  % stepsize control
  %
  if cds.h < cds.h_max & corrections==1 & it < 4
    cds.h = min(cds.h*cds.h_inc_fac, cds.h_max);
    %printconsole('<+>');
  end

  % closed curve?
  %
  if CheckClosed>0 & i>=CheckClosed & norm(x2-x0) < cds.h 
    i=i+1;
    xout(:,i) = xout(:,1);
    vout(:,i) = vout(:,1);
    hout(:,i) = hout(:,1);
    if ~isempty(fout)
      fout(:,i) = fout(:,1);
    end
    printconsole('Closed curve detected at step %d\n', i);
    break;
    %else
    %printconsole('step %d: %f\n', i, norm(x2-x0));
  end
  
  if gui_enabled
    if (mod(i,plotpoints)==0)
        sOutput.output(xout,s,hout,fout,ind);
        ind=[];
    end
    drawnow;
    doStop = sOutput.checkPauseResumeStop(s.index , s.msg ,i);
    if doStop~=0, npoints = i;end  %npoints = i ---> end of mainloop

end


    % shift x1,v1  %
  x1 = x2;
  v1 = v2;
end

printconsole('\n');
EndTime = clock;

% if closed curve, we have possibly less points
if npoints > i
  npoints = i;
  xout = xout(:,1:npoints);
  vout = vout(:,1:npoints);
  hout = hout(:,1:npoints);
%  if ~isempty(fout)
    fout = fout(:,1:npoints);
    % end
end

% Last point is also singular
% WM: Added a call to the default processor here
xout = xout(:,1:i);
vout = vout(:,1:i);
hout = hout(:,1:i);
fout = fout(:,1:i);
s.index = npoints;
s.label = '99';
s.msg   = 'This is the last point on the curve';
s.data  = [];
if gui_enabled; s.internallabel = []; end
[failed,f,s] = DefaultProcessor(xout(:,npoints),vout(:,npoints),s);
fout(:,i)= f;
ind=[ind i];
if gui_enabled
    sOutput.output(xout,sout(end),hout,fout,ind);
end
sout = [sout; s];

if WorkSpace
%   [xout,vout,sout]=feval(cds.curve, 'done',xout,vout,sout);
    feval(cds.curve_done);
end



string=sprintf('%.1f secs\n', etime(EndTime, StartTime));

printconsole('elapsed time  = %.1f secs\n', etime(EndTime, StartTime));
printconsole('npoints curve = %d\n', npoints);
%sOutput Duraction
sOutput.endRun();


%--< END OF CONTINUER >--



%-----------------------------------
%
% Command line parser
%
%-----------------------------------

function [curve, x, v, options] = ParseCommandLine(curve, x0, varargin)

% if nargin < 2 | ~isstr(curve) | isempty(x0) 
if nargin < 2 |  isempty(x0) 
  error('error in parameters, see help');
end

x = x0;
v = [];
options = contset;

if nargin > 2
  v = varargin{1};
  if nargin > 3
    options = contmrg(options,varargin{2});
  end
end
%--< END OF CMDL PARSER >--


%-------------------------------------
%
% Start point Corrector
%
%-------------------------------------

function [x,v] = CorrectStartPoint(x0, v0)
global cds

x = [];
v = [];

% no tangent vector given, cycle through base-vectors

if cds.options.TSearchOrder==1;
i = 1;
v0 = zeros(cds.ndim,1);
while isempty(x) & i<=cds.ndim
  v0(i) = 1;
   try
  DefaultProcessor(x0,v0);
  [x,v] = newtcorr(x0, v0);
   catch
   end
  v0(i) = 0; 
  i=i+1;
end
  else
i = cds.ndim;
v0 = zeros(cds.ndim,1);
while isempty(x) & i>=1
  v0(i) = 1;
   try
  DefaultProcessor(x0,v0);
  [x,v] = newtcorr(x0, v0);
   catch
   end
  v0(i) = 0; 
  i=i-1; 
end
end

%printconsole('\n');
% if ~isempty(x)
%   printconsole('Start point corrected with base vector e(%d)\n', i);
%end

%--< END OF ST PNT CORRECTOR >--



%------------------------------------------
%
%  Evaluate testfunctions
%
%------------------------------------------

function [out,failed] = EvalTestFunc(id,x,v)
global cds

if id == 0
  % WM: evaluate all testfunctions at once
  [out,failed] = feval(cds.curve_testf, 1:cds.nTest, x, v);
else
  [out,failed] = feval(cds.curve_testf, id, x, v);
end




%--< END OF TESTF EVAL >--


%------------------------------------------------------------
%
% Locate singularity xs between x1 and x2
% First locating zeros of testfunctions, then glue them
%
%------------------------------------------------------------

function [xs,vs] = LocateSingularity(si)
global cds
% zero locations of testf i is stored in testzero(:,i)

% if 1 zero/tf check if nonzero/tf _is_ nonzero at that point
% if more zero/tf then glue, nonzero/tf is kind of a problem because can always be found
idx = find( cds.S(si,:)==0 );
nzs = find( cds.S(si,:)==1 );
len = length(idx);
lnz = length(nzs);

switch len
  case 0
  % Oops, we have detected a singularity without a vanishing testfunction
  error('Internal error: trying to locate a non-detected singularity');

  case 1
  % check all nonzero/tf
  xs = cds.testzero(:,idx);
  vs = cds.testvzero(:,idx);

  otherwise
  tz = zeros(cds.ndim,len);
  vz = zeros(cds.ndim,len);
  nm = zeros(1,len);

  for i=1:len
    tz(:,i) = cds.testzero(:,idx(i));
    vz(:,i) = cds.testvzero(:,idx(i));
    nm(i) = norm(tz(:,i));
  end
 
  if max(nm)-min(nm) < 10000*cds.options.VarTolerance
    xs = mean(tz',1)';
    vs = mean(vz',1)';
  else
    xs = [];
    vs = [];
    return;
  end
end
  
if lnz == 0, return; end

% checking non zeros

DefaultProcessor(xs,vs);
tval = EvalTestFunc(nzs, xs, vs);
if any(abs(tval(nzs)) <= cds.options.TestTolerance)
  xs = [];
  vs = [];
end
%printconsole('nz of tf %d = %f\n', nzs(ti), tv);


%------------------------------------------------
%
%  Sorts [xs,vs] and sing with x1 starting point
%  xs contains x_singular, vs v_singular
%  sing contains their id's
%  Sort criterion: dist(x1,x)
%
%------------------------------------------------

function [xs,vs,sing] = xssort(x1, xs, vs, sing)
% WM: Matlab has a sort function, beter use it...
len = size(xs,2);
if len > 1
  xo = x1(:,ones(1,len));
  [dummy,i] = sortrows((xs-xo)');
  xs = xs(:,i);
  vs = vs(:,i);
  sing = sing(:,i);
end

%--< END OF XSSORT >--

%----------------------------
%
% DefaultProcessor
%
%----------------------------

function [failed,f,s] = DefaultProcessor(x,v,s)
global cds
% WM: this now actually calls the default processor,
% either with or without a singular point structure

if nargin > 2
  [failed,f,s] = feval(cds.curve_defaultprocessor, x, v, s);
else
  [failed,f] = feval(cds.curve_defaultprocessor, x, v);
end


%--< END OF DEFAULTPROCESSOR >--

%----------------------------------------------
function [x,v,i] = LocateTestFunction(id,x1,v1,x2,v2)
% default locator: bisection
global cds
%printconsole('locating tfz %d\n', id);
% WM: eliminated found variable
i = 1;
t1 = EvalTestFunc(id,x1,v1);
t2 = EvalTestFunc(id,x2,v2);
tmax = 10*max(abs(t1(id)),abs(t2(id)));
p = 1;
try
while i<=cds.options.MaxTestIters
    
  % WM: make educated guess of where the zero point might be
    if tmax < Inf
        r = abs(t1(id)/(t1(id)-t2(id)))^p;
    else
        r=0.5;
    end
    x3 = x1 + r*(x2-x1);
    v3 = v1 + r*(v2-v1);
    [x,v] = newtcorr(x3,v3);
    if isempty(x)
        x = x3;
        v = v3;
    end

  DefaultProcessor(x,v);
  tval = EvalTestFunc(id,x,v);
  dist1 = norm(x-x1);
  dist2 = norm(x-x2);
  if abs(tval(id)) > tmax
%     printconsole('testfunction behaving badly.\n');
    x = [];
    break;
  end
    
if abs(tval(id)) <= cds.options.TestTolerance & min(dist1,dist2) < 1000*cds.options.VarTolerance
    
    break;
    % elseif sign(tval(id))~=sign(cds.testvals(cds.atv,id))
    % change by WG
  elseif sign(tval(id))==sign(t2(id))
    x2 = x;
    v2 = v;
    t2 = tval;
    p = 1.02;
  else
    x1 = x;
    v1 = v;
    t1 = tval;
    p = 0.98;
  end
  i = i+1;
  x=[];
  v=[];
end
catch
    x = [];
    v = [];
    i = cds.options.MaxTestIters;
end
%--< END OF locatetestfunction>--
%---------------------------------------------
%----------------------------------------------
%
%LocateUserFunction(id,x1,v1,x2,v2)
%
%----------------------------------------------
function [x,v,i] = LocateUserFunction(userinf,id,x1,v1,x2,v2)
% default locator: bisection
global cds
%printconsole('locating tfz %d\n', id);

% WM: eliminated found variable
i = 1;
t1 = feval(cds.curve_userf, userinf, id, x1, v1);
t2 = feval(cds.curve_userf, userinf, id, x2, v2);
tmax = 10*max(abs(t1),abs(t2));
p = 1;
try
while i<=cds.options.MaxTestIters
if tmax < Inf
  % WM: make educated guess of where the zero point might be
  r = abs(t1/(t1-t2))^p;
else
    r=0.5;
end
%  r = 0.5; %  -> 'normal' way
  x3 = x1 + r*(x2-x1);
  v3 = v1 + r*(v2-v1);
  [x,v] = newtcorr(x3,v3);
  if isempty(x)
    x = x3;
    v = v3;
  end

  DefaultProcessor(x,v);
  tval=feval(cds.curve_userf,userinf,id,x,v);
  dist1 = norm(x-x1);
  dist2 = norm(x-x2);
  if abs(tval) > tmax
%     printconsole('testfunction behaving badly.\n');
    x = [];
    break;
  end
  if abs(tval) <= cds.options.TestTolerance & min(dist1,dist2) < 1000*cds.options.VarTolerance
%  if abs(tval) < cds.options.TestTolerance & dist < cds.options.VarTolerance
    break;
    % elseif sign(tval(id))~=sign(cds.testvals(cds.atv,id))
    % change by WG
elseif sign(tval)==sign(t2)
    x2 = x;
    v2 = v;
    t2 = tval;
    p = 1.02;
  else
    x1 = x;
    v1 = v;
    t1 = tval;
    p = 0.98;
  end
  i = i+1;
  x=[];  
end
catch
end
%--< END OF locateuserfunction>--

%SD:actual continuer code
