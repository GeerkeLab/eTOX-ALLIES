function displayOL(enabled,idObj) {
	if (enabled) {
		document.getElementById(idObj).style.visibility = "visible";
		document.getElementById(idObj).style.overflow = "visible";
		document.getElementById('minus'+idObj).style.display = "";
     	document.getElementById("plus"+idObj).style.display = "none"; 
  	} else {
     	document.getElementById(idObj).style.visibility = "collapse";
    	document.getElementById(idObj).style.overflow = "hidden";
     	document.getElementById("minus"+idObj).style.display = "none";
     	document.getElementById("plus"+idObj).style.display = "";
  	}
}
function showMol(smi,link) {
	smiMod=smi.replace(/=/g,'&equal&');
	smiMod=smiMod.replace(/#/g,'&hash&');
	smiMod=smiMod.replace(/\//g,'&slash&');

	linkMod=link+smiMod;
	

    var myWindow = window.open(linkMod, "MsgWindow", "width=500, height=300");
}

function changeColor(elmnt,clr) {
    elmnt.style.color = clr;
}