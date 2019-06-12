function [finalTriplets,v,testPassed,origgTime,currentTime] = buildTripletsAviodCollinearFast(pointMatchesInliers,FN)

numcam=size(pointMatchesInliers,2);
epipolsx=zeros(numcam,numcam);
epipolsy=zeros(numcam,numcam);


inliersNum= pointMatchesInliers;

adjec=zeros(numcam);
for i=1:numcam-1
    for j=i+1:numcam
        curF=FN(3*i-2:3*i,3*j-2:3*j);
        ncurF=null(curF);
        ncurFt=null(curF');
        ncurF=ncurF/ncurF(3);
        
        epipolsx(i,j)=ncurF(1);
        epipolsy(i,j)=ncurF(2);
        
        ncurFt=ncurFt/ncurFt(3);
        epipolsx(j,i)=ncurFt(1);
        epipolsy(j,i)=ncurFt(2);
        if  inliersNum(i,j)>0
            adjec(i,j)=1/inliersNum(i,j);
            adjec(j,i)=1/inliersNum(i,j);
        end
    end
end
G=graph(adjec);
center=[1296,1936]'/2;

T=minspantree(G);
tripletGraph=adjacency(T);
graphsTrips=cell(4,1);
graphsTrips{1}=graph(tripletGraph);
try
    for i=1:5
        X=setdiff(G.Edges,T.Edges);
        G=graph(X.EndNodes(:,1) , X.EndNodes(:,2),X.Weight);
        T=minspantree(G);
        tripletGraph=tripletGraph+adjacency(T);
        graphsTrips{i+1}=graph(tripletGraph);
    end
catch
    
end
[c]=extractTripletsFromViewingGraph(graph(tripletGraph));

tts=zeros(size(c,1),1);

tripletsErrors=zeros(size(c,1),2);
for i=1:size(c,1)
    e12=[epipolsx(c(i,2),c(i,1));epipolsy(c(i,2),c(i,1))];
    e13=[epipolsx(c(i,3),c(i,1));epipolsy(c(i,3),c(i,1))];
    e21=[epipolsx(c(i,1),c(i,2));epipolsy(c(i,1),c(i,2))];
    e23=[epipolsx(c(i,3),c(i,2));epipolsy(c(i,3),c(i,2))];
    e31=[epipolsx(c(i,1),c(i,3));epipolsy(c(i,1),c(i,3))];
    e32=[epipolsx(c(i,2),c(i,3));epipolsy(c(i,2),c(i,3))];
    tts(i)=getCollinearityMeasurement( e12,e13,e21,e23,e31,e32,center );
    inds=[3*c(i,1)-2:3*c(i,1) 3*c(i,2)-2:3*c(i,2) 3*c(i,3)-2:3*c(i,3) ];
    curFF=FN(inds,inds);
    curFF=normalizeF(curFF,[2000 0 1000;0 2000 1000;0 0 1]);
     [ curFF ] = normalizeForbineusNorm( curFF );
    [tripletsEr]=TripletError(curFF );
    tripletsErrors(i,1)=tripletsEr;
    tripletsErrors(i,2)=tripletsEr;
 end
c=c(tts>0.03 ,:);
tripletsErrors=tripletsErrors(tts>0.03 ,:);

tts=tts(tts>0.03,:);

adjec=zeros(size(c,1));
for i=1:size(c,1)-1
    
    for j=i+1:size(c,1)
        
        tuple1=c(i,:);
        tuple2=c(j,:);
        
        numequals=0;
        indd1=1;
        indd2=1;
        
        for k=1:4
            if tuple1(indd1)==tuple2(indd2)
                numequals=numequals+1;
                indd1=indd1+1;
                indd2=indd2+1;
            elseif  tuple1(indd1)>tuple2(indd2)
                indd2=indd2+1;
            else
                indd1=indd1+1;
            end
            if indd1>3 || indd2>3
                break;
            end
            
        end
        if numequals>=2
            adjec(i,j)=1;
            adjec(j,i)=1;
        end
      
    end
end
Gt=graph(adjec);

if mean(tts)<0.5
    tts2=(tripletsErrors(:,2).^(-1)).*((tts).^1.2);
else
    tts2=(tripletsErrors(:,2).^(-1));
end

[ firstGroup,Gt]=makeMinimalGraph([c,repmat(tts2,1,3)],1:size(c,1),Gt,numcam);

testPassed=true;

finalTriplets=c(firstGroup,:);
v = bfsearch(Gt,1,{'edgetonew'});

end