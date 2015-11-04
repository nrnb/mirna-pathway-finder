window.addEventListener('load', function() {
  function updateWidget(identifier) {
    var widget = d3.select('#wikipathways-widget').selectAll('iframe')
      .data([3])
      .attr('src', 'http://www.wikipathways.org/wpi/PathwayWidget.php?id=' + identifier + '&' + highlightString);

    widget.enter().append('iframe')
      .attr('src', 'http://www.wikipathways.org/wpi/PathwayWidget.php?id=' + identifier + '&' + highlightString)
       .attr('width', '600px')
       .attr('height', '300px')
       .attr('style', 'overflow:hidden;');

    widget.exit().remove();
  }

  d3.select('#pathway-to-mirna')
    .selectAll('tr')
    .select('#name')
    .filter(function(d) {
      return this;
    })
    .on('click', function(event) {
      var identifier = d3.select(this.parentNode).select('#identifier').text();
      console.log('identifier');
      console.log(identifier);
      updateWidget(identifier);
    });
});
