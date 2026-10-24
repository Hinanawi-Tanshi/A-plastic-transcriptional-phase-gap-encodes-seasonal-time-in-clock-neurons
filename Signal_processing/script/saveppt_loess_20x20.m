function saveppt_loess_20x20(p)
import mlreportgen.ppt.*

[~, name] = fileparts(p);
ppt = Presentation([name '.pptx']);
titleSlide = add(ppt,'Title and Content');


paraObj = Paragraph(p);
paraObj.Font = 'Arial';
paraObj.FontSize = '18pt';
replace(titleSlide,'Title',paraObj);

pictureSlide0 = add(ppt,'Blank');
plot0 = Picture('cell_labels.tif');
add(pictureSlide0,plot0);

slideAdded1 = add(ppt,'Title Slide');
replace(slideAdded1,'Title','Original Results');
pictureSlide1 = add(ppt,'Blank');
plot1 = Picture([p '\SignalFigs\IndividualSignal.png']);
add(pictureSlide1,plot1);
pictureSlide2 = add(ppt,'Blank');
plot2 = Picture([p '\SignalFigs\SignalFigures.png']);
add(pictureSlide2,plot2);

slideAdded1_raw = add(ppt,'Title Slide');
replace(slideAdded1_raw,'Title','Original Results_raw');
pictureSlide1_raw = add(ppt,'Blank');
plot1_raw = Picture([p '\SignalFigs\IndividualSignal_raw.png']);
add(pictureSlide1_raw,plot1_raw);
pictureSlide2_raw = add(ppt,'Blank');
plot2_raw = Picture([p '\SignalFigs\SignalFigures_raw.png']);
add(pictureSlide2_raw,plot2_raw);

slideAdded1_loess = add(ppt,'Title Slide');
replace(slideAdded1_loess,'Title','Original Results_loess');
pictureSlide1_loess = add(ppt,'Blank');
plot1_loess = Picture([p '\SignalFigs\IndividualSignal_loess.png']);
add(pictureSlide1_loess,plot1_loess);
pictureSlide2_loess = add(ppt,'Blank');
plot2_loess = Picture([p '\SignalFigs\SignalFigures_loess.png']);
add(pictureSlide2_loess,plot2_loess);

slideAdded2 = add(ppt,'Title Slide');
replace(slideAdded2,'Title','Max All Projection');
pictureSlide3 = add(ppt,'Blank');
plot3 = Picture([p '\SignalFigs\IndividualSignal_max-all-proj.png']);
add(pictureSlide3,plot3);
pictureSlide4 = add(ppt,'Blank');
plot4 = Picture([p '\SignalFigs\SignalFigures_max-all-proj.png']);
add(pictureSlide4,plot4);

slideAdded2_loess = add(ppt,'Title Slide');
replace(slideAdded2_loess,'Title','Max All Projection_loess');
pictureSlide3_loess = add(ppt,'Blank');
plot3_loess = Picture([p '\SignalFigs\IndividualSignal_max-all-proj_loess.png']);
add(pictureSlide3_loess,plot3_loess);
pictureSlide4_loess = add(ppt,'Blank');
plot4_loess = Picture([p '\SignalFigs\SignalFigures_max-all-proj_loess.png']);
add(pictureSlide4_loess,plot4_loess);

slideAdded3 = add(ppt,'Title Slide');
replace(slideAdded3,'Title','Mean-all-projection');
pictureSlide5 = add(ppt,'Blank');
plot5 = Picture([p '\SignalFigs\IndividualSignal_mean-all-proj.png']);
add(pictureSlide5,plot5);
pictureSlide6 = add(ppt,'Blank');
plot6 = Picture([p '\SignalFigs\SignalFigures_mean-all-proj.png']);
add(pictureSlide6,plot6);

slideAdded3_loess = add(ppt,'Title Slide');
replace(slideAdded3_loess,'Title','Mean-all-projection_loess');
pictureSlide5_loess = add(ppt,'Blank');
plot5_loess = Picture([p '\SignalFigs\IndividualSignal_mean-all-proj_loess.png']);
add(pictureSlide5_loess,plot5_loess);
pictureSlide6_loess = add(ppt,'Blank');
plot6_loess = Picture([p '\SignalFigs\SignalFigures_mean-all-proj_loess.png']);
add(pictureSlide6_loess,plot6_loess);

slideAdded4 = add(ppt,'Title Slide');
replace(slideAdded4,'Title','Max-3slice-projection');
pictureSlide7 = add(ppt,'Blank');
plot7 = Picture([p '\SignalFigs\IndividualSignal_max-3slice-proj.png']);
add(pictureSlide7,plot7);
pictureSlide8 = add(ppt,'Blank');
plot8 = Picture([p '\SignalFigs\SignalFigures_max-3slice-proj.png']);
add(pictureSlide8,plot8);

slideAdded4_loess = add(ppt,'Title Slide');
replace(slideAdded4_loess,'Title','Max-3slice-projection_loess');
pictureSlide7_loess = add(ppt,'Blank');
plot7_loess = Picture([p '\SignalFigs\IndividualSignal_max-3slice-proj_loess.png']);
add(pictureSlide7_loess,plot7_loess);
pictureSlide8_loess = add(ppt,'Blank');
plot8_loess = Picture([p '\SignalFigs\SignalFigures_max-3slice-proj_loess.png']);
add(pictureSlide8_loess,plot8_loess);

slideAdded5 = add(ppt,'Title Slide');
replace(slideAdded5,'Title','Mean-3slice-projection');
pictureSlide9 = add(ppt,'Blank');
plot9 = Picture([p '\SignalFigs\IndividualSignal_mean-3slice-proj.png']);
add(pictureSlide9,plot9);
pictureSlide10 = add(ppt,'Blank');
plot10 = Picture([p '\SignalFigs\SignalFigures_mean-3slice-proj.png']);
add(pictureSlide10,plot10);

slideAdded5_loess = add(ppt,'Title Slide');
replace(slideAdded5_loess,'Title','Mean-3slice-projection_loess');
pictureSlide9_loess = add(ppt,'Blank');
plot9_loess = Picture([p '\SignalFigs\IndividualSignal_mean-3slice-proj_loess.png']);
add(pictureSlide9_loess,plot9_loess);
pictureSlide10_loess = add(ppt,'Blank');
plot10_loess = Picture([p '\SignalFigs\SignalFigures_mean-3slice-proj_loess.png']);
add(pictureSlide10_loess,plot10_loess);

slideAdded6 = add(ppt,'Title Slide');
replace(slideAdded6,'Title','Max-3slice-projection & Max-1slice-calculation');
pictureSlide11 = add(ppt,'Blank');
plot11 = Picture([p '\SignalFigs\IndividualSignal_max-3slice-proj&max-1slice-cal.png']);
add(pictureSlide11,plot11);
pictureSlide12 = add(ppt,'Blank');
plot12 = Picture([p '\SignalFigs\SignalFigures_max-3slice-proj&max-1slice-cal.png']);
add(pictureSlide12,plot12);

slideAdded6_loess = add(ppt,'Title Slide');
replace(slideAdded6_loess,'Title','Max-3slice-projection & Max-1slice-calculation_loess');
pictureSlide11_loess = add(ppt,'Blank');
plot11_loess = Picture([p '\SignalFigs\IndividualSignal_max-3slice-proj&max-1slice-cal_loess.png']);
add(pictureSlide11_loess,plot11_loess);
pictureSlide12_loess = add(ppt,'Blank');
plot12_loess = Picture([p '\SignalFigs\SignalFigures_max-3slice-proj&max-1slice-cal_loess.png']);
add(pictureSlide12_loess,plot12_loess);

slideAdded7 = add(ppt,'Title Slide');
replace(slideAdded7,'Title','Normalized Original Results');
pictureSlide13 = add(ppt,'Blank');
plot13 = Picture([p '\SignalFigs\NormalizedIndividualSignal.png']);
add(pictureSlide13,plot13);
pictureSlide14 = add(ppt,'Blank');
plot14 = Picture([p '\SignalFigs\NormalizedSignalFigures.png']);
add(pictureSlide14,plot14);

slideAdded8 = add(ppt,'Title Slide');
replace(slideAdded8,'Title','Normalized Original Results_loess');
pictureSlide15 = add(ppt,'Blank');
plot15 = Picture([p '\SignalFigs\NormalizedIndividualSignal_loess.png']);
add(pictureSlide15,plot15);
pictureSlide16 = add(ppt,'Blank');
plot16 = Picture([p '\SignalFigs\NormalizedSignalFigures_loess.png']);
add(pictureSlide16,plot16);

slideAdded9 = add(ppt,'Title Slide');
replace(slideAdded9,'Title','Normalized Original Results: all 20x20');
pictureSlide17 = add(ppt,'Blank');
plot17 = Picture([p '\SignalFigs\NormalizedIndividualSignal_all20x20.png']);
add(pictureSlide17,plot17);

slideAdded10 = add(ppt,'Title Slide');
replace(slideAdded10,'Title','Normalized Original Results: partial 20x20');
pictureSlide18 = add(ppt,'Blank');
plot18 = Picture([p '\SignalFigs\NormalizedIndividualSignal_partial20x20.png']);
add(pictureSlide18,plot18);

close(ppt);
fprintf('Presentation saved: %s\n',fullfile(p,[name '.pptx']));
end