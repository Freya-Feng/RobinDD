function percentage_count(i,n,c)
% show percentage on srean
% i step, n total count, c name

      jd=int16(i/n*100);
      if i==1
         fprintf([c,': %3i%%'], jd)
      else      
         fprintf('\b\b\b\b')
         fprintf('%3i%%', jd)
      end
      if i==n
          fprintf('\n')
      end
    
return