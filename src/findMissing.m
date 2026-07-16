function findMissing(data,thres)
    curr = data(1);
    for i = 1:length(data)-1
        if curr ~= data(i)
            fprintf("Jump at idx="+i+"\n");
            fprintf("From: "+data(i-1)+" to "+data(i)+"\n")
            curr = data(i);
        end
        if curr==thres
            curr=0;
        else
        curr = curr+1;
        end
    end
end