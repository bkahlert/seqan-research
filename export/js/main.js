(function($){

$(function() {

$('[href]').each(function() {
  var $this = $(this);
  var href = $this.attr('href') + "";
  if(href.match('apiua://')) {
    $this.attr('href', '#' + encodeURI($this.attr('href')));
  }
});

var fadeStart=100;
var fadeUntil=300;

$(window).bind('scroll', function(){
    var offset = $(document).scrollTop()
        ,opacity=0;
    if( offset<=fadeStart ){
        opacity=1;
    }else if( offset<=fadeUntil ){
        opacity=1-offset/fadeUntil;
    }
    $('.navbar[role=navigation]').css('opacity',opacity);
});

});

})(jQuery);