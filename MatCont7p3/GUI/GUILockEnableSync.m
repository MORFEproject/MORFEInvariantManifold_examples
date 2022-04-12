classdef GUILockEnableSync < handle


   properties
       handlelist
       predicate
       listener
       listener2
   end

   methods
       function obj = GUILockEnableSync(model, handlelist, eventname, predicate)
           obj.handlelist = handlelist;
           obj.predicate = predicate;
           
           obj.listener = model.addlistener(eventname , @(src,ev) obj.sync());
           obj.listener2 = model.addlistener('shutdown' , @(src,ev) obj.destructor());
           sync(obj);
       end
       
       function sync(obj)
          enable = CLbool2text(obj.predicate());
          set(obj.handlelist, 'Enable' , enable); 
       end
       
       function destructor(obj)
          delete(obj.listener);
          delete(obj.listener2);
          delete(obj);
       end
   end
   
end 
