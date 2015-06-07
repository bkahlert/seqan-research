(function($){

$(function() {

$('[href]').each(function() {
  var $this = $(this);
  var href = $this.attr('href') + "";
  if(href.match('apiua://')) {
    $this.attr('href', '#' + encodeURI($this.attr('href')));
  }
});

});

})(jQuery);