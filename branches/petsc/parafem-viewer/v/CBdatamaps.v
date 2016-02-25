// ColorBrewer CB-SEQ-5-YlOrRd
// RGB:        | HSV:
// 189 000 038 | 
// 240 059 032 | 
// 253 141 060 | 
// 254 204 092 | 
// 255 255 178 | 

macro+Datamap CBSEQ5YlOrRd {
  float+nres ratio => (DataRange[0].UIMaxValue - DataRange[0].UIMinValue) / (dataMax - dataMin);
  DMAP.DatamapValue DatamapValue[5];
  !DatamapValue[0] { 
    v1 = 1.0; 
    v2 = 1.0; 
    v3 = 1.0; 
    v4 = 0.698;
  };
  !DatamapValue[1] { 
    v1 = 1.0; 
    v2 = 0.996; 
    v3 = 0.8; 
    v4 = 0.361;
  };
  !DatamapValue[2] { 
    v1 = 1.0; 
    v2 = 0.992; 
    v3 = 0.553; 
    v4 = 0.235;
  };
  !DatamapValue[3] { 
    v1 = 1.0; 
    v2 = 0.941; 
    v3 = 0.231; 
    v4 = 0.125;
  };
  !DatamapValue[4] { 
    v1 = 1.0; 
    v2 = 0.741; 
    v3 = 0.0; 
    v4 = 0.149;
  };
  
  /* FOUR ranges */
  DefaultDataRange DataRange[4];
  !DataRange[0] {
    size = 64;
    controlPoints => {<-.DatamapValue[0], <-.DatamapValue[1]};
    
    DataMinValue => 0 + <-.dataMin;
    DataMaxValue => <-.dataMin + ((<-.dataMax - <-.dataMin) * 0.25);
    UIMinValue = 0.0;
    UIMaxValue = 63.75;
    minActive = 1;
    maxActive = 1;
    sizeActive = 1;
    selectValues = 0;
    selectAlphaRange = 0;
    selectColorRange = 0;
  };
  !DataRange[1] {
    size = 64;
    controlPoints => {<-.DatamapValue[1], <-.DatamapValue[2]};
    
    DataMinValue => DataRange[0].DataMaxValue;
    DataMaxValue => <-.dataMin + ((<-.dataMax - <-.dataMin) * 0.5);
    UIMinValue = 63.75;
    UIMaxValue = 127.5;
    minActive = 1;
    maxActive = 1;
    sizeActive = 1;
    selectValues = 0;
    selectAlphaRange = 0;
    selectColorRange = 0;
  };
  !DataRange[2] {
    size = 64;
    controlPoints => {<-.DatamapValue[2], <-.DatamapValue[3]};
    
    DataMinValue => DataRange[1].DataMaxValue;
    DataMaxValue => <-.dataMin + ((<-.dataMax - <-.dataMin) * 0.75);
    UIMinValue = 127.5;
    UIMaxValue = 191.25;
    minActive = 1;
    maxActive = 1;
    sizeActive = 1;
    selectValues = 0;
    selectAlphaRange = 0;
    selectColorRange = 0;
  };
  !DataRange[3] {
    size = 64;
    controlPoints => {<-.DatamapValue[3], <-.DatamapValue[4]};
    
    DataMinValue => DataRange[2].DataMaxValue;
    DataMaxValue => 0 + <-.dataMax;
    UIMinValue = 191.25;
    UIMaxValue = 255.0;
    minActive = 1;
    maxActive = 1;
    sizeActive = 1;
    selectValues = 0;
    selectAlphaRange = 0;
    selectColorRange = 0;
  };

  editable = 1;
  /* allows merge operation to work properly. */
  dataMin =;
  dataMax =;
  currentColorModel = 1;
  colorModel => ColorModels.models[currentColorModel];
  ranges+IPort => DataRange;
  
  /* For compatibility, we need to have minvalue and
  maxvalue around so if people had actually modified
  them, they will be able to access them.
  */
  DMAP.DatamapValue &minvalue => DatamapValue[0];
  DMAP.DatamapValue &maxvalue => DatamapValue[4];
};

macro+Datamap CBSEQ5YlOrRdSTEP {
  float+nres ratio => (DataRange[0].UIMaxValue - DataRange[0].UIMinValue) / (dataMax - dataMin);
  DMAP.DatamapValue DatamapValue[5];
  !DatamapValue[0] { 
    v1 = 1.0; 
    v2 = 1.0; 
    v3 = 1.0; 
    v4 = 0.698;
  };
  !DatamapValue[1] { 
    v1 = 1.0; 
    v2 = 0.996; 
    v3 = 0.8; 
    v4 = 0.361;
  };
  !DatamapValue[2] { 
    v1 = 1.0; 
    v2 = 0.992; 
    v3 = 0.553; 
    v4 = 0.235;
  };
  !DatamapValue[3] { 
    v1 = 1.0; 
    v2 = 0.941; 
    v3 = 0.231; 
    v4 = 0.125;
  };
  !DatamapValue[4] { 
    v1 = 1.0; 
    v2 = 0.741; 
    v3 = 0.0; 
    v4 = 0.149;
  };
  
  /* FOUR ranges */
  DataRange DataRange[1];
  !DataRange[0] {
    size = 5;
    controlPoints => {<-.DatamapValue[0], <-.DatamapValue[1], 
    <-.DatamapValue[2],<-.DatamapValue[3],<-.DatamapValue[4]};
    
    DataMinValue => 0 + <-.dataMin;
    DataMaxValue => 0 + <-.dataMax;
    UIMinValue = 0.0;
    UIMaxValue = 255.0;
    minActive = 1;
    maxActive = 1;
    sizeActive = 1;
    selectValues = 0;
    selectAlphaRange = 1;
    selectColorRange = 1;
  };

  editable = 1;
  /* allows merge operation to work properly. */
  dataMin =;
  dataMax =;
  currentColorModel = 1;
  colorModel => ColorModels.models[currentColorModel];
  ranges+IPort => DataRange;
  
  /* For compatibility, we need to have minvalue and
  maxvalue around so if people had actually modified
  them, they will be able to access them.
  */
  DMAP.DatamapValue &minvalue => DatamapValue[0];
  DMAP.DatamapValue &maxvalue => DatamapValue[4];
};
