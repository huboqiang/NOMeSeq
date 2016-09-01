def AddHead(fid,GraphID):
    ''' AddHead     Add the head content of an XGMML file to a file handle. '''
    
    import time;
    fid.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
    fid.write('<graph id="423" label="436" directed="1" cy:documentVersion="3.0" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">\n')
    fid.write('\t<att name="networkMetadata">\n')
    fid.write('\t\t<rdf:RDF>\n')
    fid.write('\t\t\t<rdf:Description rdf:about="http://www.cytoscape.org/">\n')
    fid.write('\t\t\t\t<dc:type>Protein-Protein Interaction</dc:type>\n')
    fid.write('\t\t\t\t<dc:description>N/A</dc:description>\n')
    fid.write('\t\t\t\t<dc:identifier>N/A</dc:identifier>\n')
    fid.write('\t\t\t\t<dc:date>%d-%d-%d %d:%d:%d</dc:date>\n' % time.localtime()[:6])
    fid.write('\t\t\t\t<dc:title>436</dc:title>\n')
    fid.write('\t\t\t\t<dc:source>http://www.cytoscape.org/</dc:source>\n')
    fid.write('\t\t\t\t<dc:format>Cytoscape-XGMML</dc:format>\n')
    fid.write('\t\t\t</rdf:Description>\n')
    fid.write('\t\t</rdf:RDF>\n')
    fid.write('\t</att>\n')
    fid.write('\t<att name="shared name" value="%s" type="string"/>\n' % GraphID)
    fid.write('\t<att name="selected" value="1" type="boolean"/>\n')
    fid.write('\t<att name="name" value="%s" type="string"/>\n' % GraphID)
    fid.write('\t<att name="__Annotations" type="list">\n')
    fid.write('\t</att>\n')
    fid.write('\t<graphics>\n')
    fid.write('\t\t<att name="NETWORK_WIDTH" value="1000.0" type="string"/>\n')     # Network_Width to be specified. 
    fid.write('\t\t<att name="NETWORK_NODE_SELECTION" value="true" type="string"/>\n')
    fid.write('\t\t<att name="NETWORK_CENTER_X_LOCATION" value="500.0" type="string"/>\n')
    fid.write('\t\t<att name="NETWORK_TITLE" value="" type="string"/>\n')
    fid.write('\t\t<att name="NETWORK_EDGE_SELECTION" value="true" type="string"/>\n')
    fid.write('\t\t<att name="NETWORK_DEPTH" value="0.0" type="string"/>\n')
    fid.write('\t\t<att name="NETWORK_CENTER_Z_LOCATION" value="0.0" type="string"/>\n')
    fid.write('\t\t<att name="NETWORK_SCALE_FACTOR" value="1.0" type="string"/>\n')
    fid.write('\t\t<att name="NETWORK_BACKGROUND_PAINT" value="#ffffff" type="string"/>\n')     # Background_Color to be specified. 
    fid.write('\t\t<att name="NETWORK_HEIGHT" value="500.0" type="string"/>\n')     # Network_Height to be specified. 
    fid.write('\t\t<att name="NETWORK_CENTER_Y_LOCATION" value="0.0" type="string"/>\n')
    fid.write('\t</graphics>\n')
    
    return True

def AddNode(fid,Label,ID,fill='#ff0000',shape='ELLIPSE',x=0.0,y=0.0,w=35.0,h=35.0,FontSize=15,LabelColor='#000000',Transparency=255,LabelPosition='C,C,c,0,0'):
    ''' AddNode     Add the content of a new node to an XGMML file.
    
    fid             The file handle of the XGMML file.
    Label           The label of the node.
    ID              ID of the node. 
    fill            The color of the node.
    shape           The shape of the node.  {'ELLIPSE','ROUND_RECTANGLE','TRIANGLE','DIAMOND',
                                             'PARALLELOGRAM','HEXAGON','RECTANGLE','OCTAGON','V'}
    x & y           The coordinates of the node.
    w & h           The width and height of the node. 
    FontSize        The font size of the label.
    LabelColor      The color of the label.
    Transparency    The transparency of the node. 
    LabelPosition   The label position. 
    '''
    
    fid.write('\t<node id="%s" label="%s">\n' % (ID,Label) )
    fid.write('\t\t<att name="shared name" value="%s" type="string"/>\n' % Label)
    fid.write('\t\t<att name="selected" value="0" type="boolean"/>\n')
    fid.write('\t\t<att name="name" value="%s" type="string"/>\n' % Label)
    fid.write('\t\t<graphics fill="%s" x="%f" type="%s" z="0.0" outline="#333333" width="3.0" y="%f" h="%f" w="%f">\n' % (fill,x,shape,y,h,w))
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_4" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_3" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_TRANSPARENCY" value="%d" type="string"/>\n' % Transparency)
    fid.write('\t\t\t<att name="NODE_LABEL_TRANSPARENCY" value="255" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_6" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_BORDER_TRANSPARENCY" value="255" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_4" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_2" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_7" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_5" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_5" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_6" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_9" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_7" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_LABEL" value="%s" type="string"/>\n' % Label)
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_4" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_5" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_BORDER_STROKE" value="SOLID" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_8" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_VISIBLE" value="true" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_LABEL_FONT_FACE" value="Dialog.plain,plain,12" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_DEPTH" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_2" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_8" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_4" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_LABEL_WIDTH" value="200.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_LABEL_FONT_SIZE" value="%d" type="string"/>\n' % FontSize)
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_3" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_1" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_LABEL_COLOR" value="%s" type="string"/>\n' % LabelColor)
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_1" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_TOOLTIP" value="" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_SELECTED_PAINT" value="#ffff00" type="string"/>\n')    # Seleted_Paint to be specified. 
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_9" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_2" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_LABEL_POSITION" value="%s" type="string"/>\n' % LabelPosition)
    fid.write('\t\t\t<att name="NODE_SELECTED" value="false" type="string"/>\n')    # is_Selected to be specified. 
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_7" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_9" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_8" value="org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ]," type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_SIZE_3" value="0.0" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_NESTED_NETWORK_IMAGE_VISIBLE" value="true" type="string"/>\n')
    fid.write('\t\t\t<att name="NODE_CUSTOMGRAPHICS_POSITION_1" value="C,C,c,0.00,0.00" type="string"/>\n')
    fid.write('\t\t</graphics>\n')
    fid.write('\t</node>\n')
    
    return True

def AddEdge(fid,source,target,ID,Label='',directed='1',SourceArrowShape='ARROW',TargetArrowShape='ARROW',LineType='SOLID',LineWidth=2.0,LineColor='#333333'):
    ''' AddEdge     Add the content of a new edge to an XGMML file.
        
        fid                 The file handle of the XGMML file.
        source              The start node of the edge.
        target              The end node of the edge.
        ID                  The ID of the edge, which must be special.
        Label               The label of the edge, which will be present on the graph but not necessarily special.
        directed            Whether the edge is directed.
        SourceArrowShape    The shape of the source arrow. {'ARROW','DELTA','CIRCLE','HALF_BOTTOM','DIAMOND','HALF_TOP','NONE','T'}
        TargetArrowShape    The shape of the target arrow. {'ARROW','DELTA','CIRCLE','HALF_BOTTOM','DIAMOND','HALF_TOP','NONE','T'}
        LineType            The presented type of the line. {'SOLID','DASHDOT','ZIGZAG','FORWARD_SLASH','SINEWAVE','SEPARATE_ARROW',
                                                             'BACKWARD_SLASH','EQUAL_DASH','PARALLEL_LINES','DASH','VERTICAL_SLASH',
                                                             'CONTIGUOUS_ARROW','DOTS'}
        LineWidth           The width of the line.
        LineColor           The color of the line. '''
        
    fid.write('\t<edge id="%s" label="%s (interaction) %s" source="%s" target="%s" cy:directed="%s">\n' % (ID,source,target,source,target,directed))
    fid.write('\t\t<att name="interaction" value="interaction" type="string"/>\n')
    fid.write('\t\t<att name="shared name" value="%s (interaction) %s" type="string"/>\n' % (source,target))
    fid.write('\t\t<att name="selected" value="0" type="boolean"/>\n')
    fid.write('\t\t<att name="name" value="%s (interaction) %s" type="string"/>\n' % (source,target))
    fid.write('\t\t<att name="shared interaction" value="interaction" type="string"/>\n')
    fid.write('\t\t<graphics width="%f" fill="%s">\n' % (LineWidth,LineColor))
    fid.write('\t\t\t<att name="EDGE_TRANSPARENCY" value="255" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_SOURCE_ARROW_SHAPE" value="%s" type="string"/>\n' % SourceArrowShape)
    fid.write('\t\t\t<att name="EDGE_SELECTED" value="false" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_TOOLTIP" value="" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_TARGET_ARROW_UNSELECTED_PAINT" value="#000000" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_LABEL_FONT_FACE" value="Dialog.plain,plain,10" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_TARGET_ARROW_SELECTED_PAINT" value="#ffff00" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_CURVED" value="true" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_LABEL_FONT_SIZE" value="10" type="string"/>\n')  # Label_Font_Size to be specified. 
    fid.write('\t\t\t<att name="EDGE_LINE_TYPE" value="%s" type="string"/>\n' % LineType)
    fid.write('\t\t\t<att name="EDGE_VISIBLE" value="true" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_SOURCE_ARROW_UNSELECTED_PAINT" value="#000000" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_STROKE_SELECTED_PAINT" value="#ff0000" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_SOURCE_ARROW_SELECTED_PAINT" value="#ffff00" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_LABEL" value="%s" type="string"/>\n' % Label)
    fid.write('\t\t\t<att name="EDGE_LABEL_TRANSPARENCY" value="255" type="string"/>\n')
    fid.write('\t\t\t<att name="EDGE_TARGET_ARROW_SHAPE" value="%s" type="string"/>\n' % TargetArrowShape)
    fid.write('\t\t\t<att name="EDGE_LABEL_COLOR" value="#000000" type="string"/>\n') # Label color to be specified. 
    fid.write('\t\t\t<att name="EDGE_BEND" value="" type="string"/>\n')
    fid.write('\t\t</graphics>\n')
    fid.write('\t</edge>\n')
    
    return True

