
	function enableItem (item,disable) {
		document.getElementById(item).disabled=disable;
	}

	function checkItem (item,disable) {
		document.getElementById(item).checked=disable;
	}


	function enableClass (className,active) {
		var itemsClass= document.getElementsByClassName(className);

		console.log(itemsClass.length);
		for( var i = 0; i<itemsClass.length; i++) {
			console.log(itemsClass[i].id);
			itemsClass[i].disabled=(! active);
		}	
	}


	function checkClass (className,active) {
		var itemsClass= document.getElementsByClassName(className);

		for( var i = 0; i<itemsClass.length; i++) {
			itemsClass[i].checked=active;
		}	
	}



	function enableVers (versions) {
		listVers=document.getElementsByClassName('modelVer');
		for (var i=0; i < listVers.length; i++) {
			listVers[i].hidden=false;
		}
		for (i in versions) {
			document.getElementById('ver'+i).hidden=false;
		}
	}
	
	function isChecked (itemId) {   	
    	return document.getElementById(itemId).checked;
	}

	function addProtInfo (ori) {
		console.log("start");
    	original = document.getElementById(ori);
    	var clone = original.cloneNode(true); // "deep" clone
    	clone.id = ori + ++nProts;
    	var children=clone.getElementsByTagName('*');
    	// or clone.id = ""; if the divs don't need an ID
    	var text = document.createElement('p');
    	text.innerHTML="<p style='font-size:120%;'><b>Protein Conformation " + nProts +"<b></p>";
    	var X = "Xcoord"+nProts;
    	var Y = "Ycoord"+nProts;
    	var Z = "Zcoord"+nProts;
    	var guessC = "guessCoords" + nProts;

    	for( var i = 0; i<children.length; i++) {
 		    if (i == 0) {
			children[i].parentNode.insertBefore(text, children[i]);
 		    }
 		    if (children[i].id=='guessCoords') {
 		    	console.log(children[i].id);
 		    	console.log(children[i].onchange);
 		    	children[i].removeAttribute("onchange");
				console.log(children[i].onchange);
				newFunct="var act=isChecked('" + guessC + "'); enableItem('"+ X +
				"',act);enableItem('"+ Y +"',act);enableItem('" + Z + "',act);";
				children[i].setAttribute("onchange",newFunct);
				console.log(children[i].onchange);
 		    }
 		    children[i].id=children[i].id + nProts;
 		    children[i].setAttribute("name", (children[i].name + nProts));
		}
    	
    	original.parentNode.appendChild(clone);
    	
    	document.getElementById('remProtBut').disabled=false;
	}
	
	function removeProtInfo (divID) {
		var divIDn=divID+nProts;
		var elem = document.getElementById(divIDn);
		console.log(divIDn);
		
		if (nProts >1) {
			--nProts;
			elem.parentNode.removeChild(elem);
			console.log(nProts);
		}
		
		if (nProts <2) {
			document.getElementById('remProtBut').disabled=true;
		}
		
	}