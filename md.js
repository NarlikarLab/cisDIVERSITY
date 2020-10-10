var linwidth=3, size=1,noOfmotifs,noOfseqs,scrnwidth;
var cellpad=2,sf,seqsUptoNow,highlightedMotif=-1;
scrnwidth = screen.width;
noOfmotifs = wholedata["motifs"];
noOfseqs = wholedata["totalseqs"];

function addContent(tbody,moduleNo,currseqno,totseqs,seqTags,canvascol,winheight,scrollHeight){
	var seq,seqid;
	var cv,ctx ;
	var i,motifid,startpos,width,yval;
	var tr,td,node;
    for(i = currseqno; (i<(currseqno+100) && (i<totseqs)); i++){
    	seq = seqTags[i];
		tr = tbody.insertRow();
		tr.id = "Row"+moduleNo+"_"+i;
      td = tr.insertCell();
      if (seqdetails[seq][0] == "undefined"){
     	  seqid = i;
      }
      else {
      	seqid = seqdetails[seq][0];
      }
      node = document.createTextNode(seqid+'.');
      node.fontSize = 2;
	   td.appendChild(node);
					
		td = tr.insertCell();
		node = document.createTextNode(seq);
		node.fontSize = 2;
		td.appendChild(node);
			
		td = tr.insertCell();
		cv = document.createElement("canvas");
		cv.width = canvascol;
		cv.height = Math.ceil(2.604*(screen.height)/100);
		cv.id = "Canvas"+moduleNo+"_"+i;
		td.appendChild(cv);
		yval = (cv.height)/2;
		ctx = cv.getContext('2d');				
				
		ctx.scale(sf,1);
		ctx.beginPath();
		ctx.moveTo(0, yval);
		
		ctx.lineTo(seqdetails[seq][1], yval);
		ctx.closePath();
		ctx.lineWidth =linwidth;
		ctx.strokeStyle = "#FAEBD7";
		ctx.closePath();
		ctx.stroke();				
		
		var motifsThisSeq = allmodules[moduleNo][seq];
		var mots = Object.keys(motifsThisSeq);
		if (mots.length == 0){
			continue;
		}
			
		var ele2 = document.getElementById("Canvas"+moduleNo+"_"+i);		
		if (ele2==null){
			document.write("This canvas element doesn't exist! "+"Canvas"+moduleNo+"_"+i);
		}			
		var rect = ele2.getBoundingClientRect();
		var motnoblk;
		
		for (var k=0;k<(mots.length);k++)
		{
			motifid = mots[k];
			width = motifwidths[motifid];
			startpos = motifsThisSeq[motifid]["pos"];	
			
			var direc = motifsThisSeq[motifid]["dir"];
			if (highlightedMotif > -1){
				if (k==highlightedMotif){
					if (moduleMotifColors[moduleNo][motifid]==="#FF0000")
						drawMotif(ctx,startpos,width,"#8B0000",yval,direc);
					else
						drawMotif(ctx,startpos,width,"#00FF00",yval,direc);
							
				}
				else drawMotif(ctx,startpos,width,moduleMotifColors[moduleNo][motifid-1],yval,direc);
			}	
			else{				
				drawMotif(ctx,startpos,width,moduleMotifColors[moduleNo][motifid-1],yval,direc);
			}
			var xcoord,ycoord,block,make_popup;
			//document.write("Rectangle: "+ rect.left+" "+rect.top);
			xcoord = rect.left + ((startpos)*sf);
			ycoord = scrollHeight+rect.top + yval-cellpad;
			//document.write("Coordinates: "+xcoord+" "+ycoord);
			block = drawDiv(width,xcoord,ycoord);
			document.body.appendChild(block);
			/*
			motnoblk = drawDivMod(15,12,(xcoord+width+5),ycoord-10);
			motnoblk.style.textAlign = "center";		
			motnoblk.style.fontSize = "12px";
			motnoblk.appendChild(document.createTextNode(""+motifid));
			document.body.appendChild(motnoblk);
			*/		
			
			make_popup = create_popup_block(xcoord,ycoord,motifid,i,direc);  
			block.addEventListener("mouseover",make_popup,false);
			block.addEventListener("mouseout",make_popup,false); 
		}
   }
   if (currseqno >= totseqs){
		window.removeEventListener('scroll',addContent,false);
		return;   
   }
   return i; 	
}

function tableCreate(moduleId,moduleNo){
		function add_to_header(row,tname,spacewidth){
			var rh,th,title;
			rh = row.insertCell();
			th = document.createElement("th");
			th.style.width=spacewidth;
			title = document.createTextNode(""+tname);
			th.appendChild(title);
			th.style.fontFamily="arial";
			rh.appendChild(th);
			row.appendChild(rh);
		}
		var element1 = document.getElementById(moduleId);
		if (element1 == null) {
			document.write("Cannot fetch");
		}
    	var tbl  = document.createElement("table");
		var rhead,thead,tbody;
		var snowidth,seqtagwidth,canvascol,currseqno,winheight,currheight;
		currseqno = 0;
		seqUptoNow = currseqno;
		element1.appendChild(tbl);	
		tbl.className = "motifs_tbl";
		tbl.id = "motifs_table_"+moduleNo;
		tbl.style.textAlign="left";
		tbl.cellPadding = cellpad+"px";
  		thead = document.createElement("thead");
      	tbl.appendChild(thead);
      	tbody = document.createElement("tbody");
      	tbl.appendChild(tbody);
      	tbl.style.width  = '100%';
		rhead = thead.insertRow();
		rhead.id ='RowHead';
		
		snowidth = Math.ceil((0.732064*scrnwidth)/100);	
		seqtagwidth = Math.floor((18.3016*scrnwidth)/100);
		canvascol = Math.ceil((73.2064*scrnwidth)/100);

		add_to_header(rhead,"No.",snowidth+"px");
		add_to_header(rhead,"Sequence Tags",seqtagwidth+"px");
		add_to_header(rhead,"Motifs on Sequences",canvascol+"px");
			
		sf = ((canvascol-4)*0.8)/maxlength;		
		var seqsThisModule = allmodules[moduleNo];
		var seqTags = Object.keys(seqsThisModule), totseqs, dsh;
		totseqs = seqTags.length;
		winheight = window.innerHeight;
		dsh = 0
		var count=0;
		currseqno = addContent(tbody,moduleNo,currseqno,totseqs,seqTags,canvascol,winheight,0);
		seqsUptoNow = currseqno;	
		window.addEventListener('scroll',function () {
			if(winheight + window.scrollY === document.documentElement.scrollHeight){	
				count=count+1;
				dsh = window.scrollY;
				currseqno=addContent(tbody,moduleNo,currseqno,totseqs,seqTags,canvascol,winheight,dsh);
				seqsUptoNow=currseqno;				
			}
		});	 
}
function drawDiv(motwidth,xcoord,ycoord){
		var block = document.createElement("div");
		var scalefactor=sf;
		block.className = "motif_block";
		block.style.position = "absolute";
		//block.style.cssFloat="left";
		block.style.left = xcoord+'px';
		block.style.top = ycoord+'px';
		block.style.width = ((motwidth)*scalefactor)+'px';
		block.style.height = (linwidth)+'px';
		//block.style.backgroundColor="black";		
		return block;
}

function create_popup_block(xcoord,ycoord,motifid,seqno,direc)
{
	return function(e){
	var popup,thisblockId,scalefactor;
	if (!e) var e= window.event;
	popup = create_popup_block.popup;
	thisblockId = "mot_"+seqno+"_"+motifid;
	if (e.type === "mouseover")
	{
		popup = document.getElementById(thisblockId);
		if (popup) {return};
		popup = getImage(motifid,direc);
		popup.id = thisblockId;
		scalefactor=sf;
		position_popup(motifid,xcoord,ycoord,scalefactor,seqno,popup);
		document.body.appendChild(popup);
	}
	else if (e.type === "mouseout")
	{
		popup = document.getElementById(thisblockId);
		if (popup){
			popup.parentNode.removeChild(popup);
			create_popup_block.popup = null;
		}
	}
	};
}

function getImage(motifid,direc){
		var motiffile,width,imgwidth,imgheight,img; 
		const zeroPad = (num, places) => String(num).padStart(places, '0');
		if (direc == "R"){
			motiffile = "sites_"+zeroPad(motifid,2)+".png";
		}
		else {
			motiffile = "revsites_"+zeroPad(motifid,2)+".png";
		}
		
		width = motifwidths[motifid];
		//imgwidth = 12 * width;
		//imgheight = 3.5 * width;
		imgwidth = logoSizes[motifid][0]*0.5;
		imgheight = logoSizes[motifid][1]*0.5;
		img = document.createElement("img");
		img.src = motiffile;
		img.width = imgwidth;
		img.height = imgheight;
		return img;
}

function position_popup(motifid,xcoord,ycoord,scalefactor,seqno,popup){
		var mh,mw,mx,my,ph,pw,p_x,p_y,space,page_w;
		var tmplnode = popup;
		space = 5;
		mx = xcoord;
		my = ycoord; 
		mh = linwidth*scalefactor;
		mw = motifwidths[motifid];
	
		pw = tmplnode.width;
		ph = tmplnode.height;
	
		if (window.innerWidth) {
   	 	page_w = window.innerWidth;
  		} else if (document.body) {
 	   	page_w = document.body.clientWidth;
		}
		p_x = mx + (mw/2) - (pw/2);
		p_y = my + mh + space;
		if ((page_w !=null) && (p_x > page_w-space)){
			p_x = page_w - space;
		}
		if (p_y < space){
			p_y = space;
		}
		if ((p_y+motifwidths[motifid]*sf)>document.documentElement.scrollHeight){
			p_y = ycoord - (motifwidths[motifid]*sf);
		}

		tmplnode.style.position = "absolute";
		tmplnode.style.left = p_x+'px';
		tmplnode.style.top = p_y+'px';
}
function drawMotif(ctx,startpos,width,motifcolor,yval,direc)
{
	var startX, startY;
	ctx.beginPath();
	if (direc=="L"){
		ctx.moveTo(startpos+(0.07*size),yval);
		ctx.lineTo((startpos+width),yval);
	}
	else {
		ctx.moveTo(startpos,yval);
		ctx.lineTo((startpos+width-(0.07*size)),yval);
	}	
	ctx.closePath();
	ctx.lineWidth=linwidth;
	ctx.strokeStyle=motifcolor;
	ctx.stroke();
        
	if (direc == "L"){
		startX=startpos;
		startY=yval;
		drawArrowLeft(ctx,startX+1,startY,motifcolor);
	}
	else if (direc == "R") {
		startX=startpos+width;
		startY =yval;
		drawArrowRight(ctx,startX,startY,motifcolor);
	}
}

function drawArrowRight(ctx,startX,startY,color)
{
	var arrowX,arrowTopY, arrowBottomY;
	arrowX= startX-0.07*size;
	arrowTopY=startY-(0.07*size);
	arrowBottomY=startY+(0.07*size);
	ctx.beginPath();
	ctx.moveTo(arrowX,arrowTopY);
	ctx.lineTo(startX,startY);
	ctx.lineTo(arrowX, arrowBottomY); 
	ctx.closePath();
	ctx.strokeStyle=color;
	ctx.lineWidth=linwidth;
	ctx.stroke();
}

function drawArrowLeft(ctx,startX,startY,color)
{
	var arrowX,arrowTopY, arrowBottomY;
	arrowX= startX+0.07*size;
	arrowTopY=startY-(0.07*size);
	arrowBottomY=startY+(0.07*size);
	ctx.beginPath();
	ctx.moveTo(arrowX,arrowTopY);
	ctx.lineTo(startX,startY);
	ctx.lineTo(arrowX, arrowBottomY); 
	ctx.closePath();
	ctx.strokeStyle=color;
	ctx.lineWidth=linwidth;
	ctx.stroke();
}

function buildModule(moduleNo){	
	var moddiv,div2, div3,tempd,seqmotifDiv,motsize;
	var rect,motpercent,para,noOfSeqs,pDiv;
	var noOfmotifs = wholedata["motifs"];
	
	//document.write("screen height and width: "+screen.height+ " "+screen.width);
	noOfSeqs = modseqs[moduleNo];
	motsize = (8*scrnwidth/10)/noOfmotifs;
	pDiv = document.getElementById("Module_"+moduleNo);
	moddiv = document.createElement("div");
	moddiv.id = "ModuleSummary_"+moduleNo;
	moddiv.style.backgroundColor="#f5ebf9";
	moddiv.style.borderRadius="5px";	
	moddiv.style.height="110px";
	moddiv.style.width= (9*(scrnwidth)/10)+"px";							
	moddiv.style.color = "black";		
	pDiv.appendChild(moddiv);
		
	div2 = document.createElement("div");
	moddiv.appendChild(div2);
	div2.className = "dropDown";
	div2.style.padding = "5px 2px";
	div2.innerHTML = "<font face=\"arial\">Module "+moduleNo+" ("+(noOfSeqs)+" sequences)"+"</font>";
	
	div3 = document.createElement("div");
	div3.className="arrow-right";
	div2.style.fontWeight = "bold";	
	div3.id="arrow"+moduleNo;
	div2.appendChild(div3);
	
	para = document.createElement("p");
	para.id = "para"+moduleNo;
	div2.appendChild(para);			
	pDiv.appendChild(document.createElement('br'));
		
	motpercent = moduleMotifPercent[moduleNo];
	for (var j=0;j<motpercent.length;j++){
		var cv,ctx,cvtext,ctxtext,yheight;
		cv = document.createElement("canvas");
		cv.id = "cv_"+moduleNo+"_"+(j+1);
		tempd = document.getElementById("para"+moduleNo);
		cv.width = motsize;
		cv.height = 60;
		ctx = cv.getContext('2d');
		tempd.appendChild(cv);
		ctx.fillStyle = moduleMotifColors[moduleNo][j];
		yheight = (motpercent[j]) *50;
		ctx.fillRect(15,(50-yheight-5),(0.35*motsize),yheight);
		var rect = cv.getBoundingClientRect();
		var xcoord, ycoord, flag;
		xcoord = rect.left+15;
		ycoord = rect.top+(50-yheight-5);
		var block = drawDivMod((0.35*motsize),yheight,xcoord,ycoord);
		block.id = "block_"+(j+1);
		document.body.appendChild(block);
		var nextblk,nexty;
		nexty = 15;
		nextblk = drawDivMod((0.35*motsize),nexty,xcoord,(ycoord+yheight));
		nextblk.style.textAlign = "center";		
		nextblk.appendChild(document.createTextNode(""+(j+1)));
		document.body.appendChild(nextblk);
		
		var make_popup=create_popup_blockMod((0.35*motsize),yheight,xcoord,ycoord,moduleNo,j+1);
		var highlight=highlightmotif(moduleNo,block.id,j+1);	
		block.addEventListener("mouseover",make_popup,false);
		block.addEventListener("mouseout",make_popup,false);
		block.addEventListener("click",highlight,false); 
      
	}
}

function drawDivMod(width,height,xcoord,ycoord,divcolor){
		var block = document.createElement("div");
		block.className = "percentblock";
		block.style.position = "absolute";
		block.style.left = xcoord+'px';
		block.style.top = ycoord+'px';
		block.style.width = width+'px';
		block.style.height = height+'px';
		if (divcolor!="") block.style.backgroundColor=divcolor;		
		return block;
	}
function create_popup_blockMod(width,height,xcoord,ycoord,moduleId,motifId)
{
	return function(e){
		var popup,thisblockId,defaultdirec;
		if (!e) var e= window.event;
		popup = create_popup_blockMod.popup;
		defaultdirec = "R"
		thisblockId = "pblock_"+moduleId+"_"+motifId;
		if (e.type === "mouseover")
		{
			popup = document.getElementById(thisblockId);
			if (popup) {return; }
			popup = getImage(motifId,defaultdirec);
			popup.id = thisblockId;
			position_popupMod(width,height,xcoord,ycoord,moduleId,popup);
			document.body.appendChild(popup);
		}
		else if (e.type === "mouseout")
		{
			popup = document.getElementById(thisblockId);
			if (popup){
				popup.parentNode.removeChild(popup);
				create_popup_blockMod.popup = null;
			}
		}
		};
}

function hexToRgb(hex) {
  var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  return result ? {
    r: parseInt(result[1], 16),
    g: parseInt(result[2], 16),
    b: parseInt(result[3], 16)
  } : null;
}

function checkifRed(col) {
	var rgbformat = hexToRgb(col);
	if (rgbformat.r == 255 && rgbformat.g == 0 && rgbformat.b == 0) return true;
	else return false;
}

function highlightmotif(moduleNo,blockid,motifId){
	return function (e) {
			if(!e) var e = window.event;
			
			var col,blk = document.getElementById("block_"+motifId);
			if (blk.className ==="percentblock"){				
				blk.className = "mycolorblock";
				//blk.style.backgroundColor = "#FFB6C1";
				highlightedMotif = motifId;
				if (checkifRed(moduleMotifColors[moduleNo][motifId-1])) {
					blk.style.backgroundColor = "#8B0000";
					col = "#8B0000";
				}
				else {
					blk.style.backgroundColor = "#00FF00";
					col = "#00FF00";
				}
				changeMotifColor(moduleNo,motifId,col);	
			} else if (blk.className==="mycolorblock") {
				blk.className = "percentblock";		
				blk.style.backgroundColor =	"transparent";
				highlightedMotif = -1;
				changeMotifColor(moduleNo,motifId,moduleMotifColors[moduleNo][motifId-1]);
			}		
	};	
}


// change the color of the motif with Id=x to color=y on all the sequences upto now
function changeMotifColor(moduleNo,thismotif,toColor) {
	var i,canvasele,ctx,yval;
	var k,seqTags,seq,rect,mots,motifid,width,direc,startpos;
	for (i=0;i<seqsUptoNow;i++){
		canvasele = document.getElementById("Canvas"+moduleNo+"_"+i);
		if (canvasele == null){ 
			document.write("This canvas element doesn't exist! "+"Canvas"+moduleNo+"_"+i);
		}
		ctx = canvasele.getContext('2d');
		yval = (canvasele.height)/2;
		seqTags = Object.keys(allmodules[moduleNo]);
		seq = seqTags[i];
		mots = Object.keys(allmodules[moduleNo][seq]);
		rect = canvasele.getBoundingClientRect();
		for (k=0;k<(mots.length);k++)
		{
			motifid = mots[k];
			width = motifwidths[motifid];
			startpos = allmodules[moduleNo][seq][motifid]["pos"];
			direc = allmodules[moduleNo][seq][motifid]["dir"];
			if (motifid==thismotif){
				drawMotif(ctx,startpos,width,toColor,yval,direc)	
			}
			//else drawMotif(ctx,startpos,width,motifColors[motifid],yval,direc);
		}
	}

}

function position_popupMod(width,height,xcoord,ycoord,moduleId,popup){
		var mh,mw,mx,my,ph,pw,p_x,p_y,space,page_w;
		var tmplnode = popup;
		space = 5;
		mx = xcoord;
		my = ycoord; 
		mh = height;
		mw = width;
	
		pw = tmplnode.width;
		ph = tmplnode.height;
	
		if (window.innerWidth) {
   	 	page_w = window.innerWidth;
  		} else if (document.body) {
 	   	page_w = document.body.clientWidth;
		}
		p_x = mx + (mw/2) - (pw/2);
		p_y = my + mh + space;
		if ((page_w !=null) && (p_x > page_w-space)){
			p_x = page_w - space;
		}
		if (p_x<space) p_x = space
		
		if (p_y < space){
			p_y = space;
		}

		tmplnode.style.position = "absolute";
		tmplnode.style.left = p_x+'px';
		tmplnode.style.top = p_y+'px';
	}
