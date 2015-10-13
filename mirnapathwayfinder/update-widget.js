function updateWidget(identifier, displayNames) {
  var highlightString = displayNames
    .map(function(displayName) {
      return 'label[]=' + encodeURIComponent(displayName);
    })
    .join('&') + '&colors=red';

  // Update
  var widget = d3.select('body').select('#wikipathways-widget').selectAll('iframe')
    .data([3])
    .attr('src', 'http://www.wikipathways.org/wpi/PathwayWidget.php?id=' + identifier + '&' + highlightString);

  // Enter
  widget.enter().append('iframe')
    .attr('src', 'http://www.wikipathways.org/wpi/PathwayWidget.php?id=' + identifier + '&' + highlightString)
     .attr('width', '600px')
     .attr('height', '300px')
     .style('overflow:hidden;');

  // Exit
  widget.exit().remove();
}

d3.select('#pathway-to-mirna')
  .selectAll('tr')
  .select('td')
  .filter(function(d) {
    return this;
  })
  .on('click', function(event) {
    var identifier = this.textContent;
    var myo = d3.select('#pathway-to-mirna')
      .selectAll('tr')
      .select('td')
      .filter(function(d) {
        return this;
      })[0][0];
    var displayNames = d3.select(myo.parentNode).selectAll('td')
      .map(function(els) {
        return els.map(function(el) {
          return el.textContent;
        })
        .filter(function(text) {
          return !!text && text.match(/^hsa/);
        })
        .map(function(text) {
          return text.split(',')
          .filter(function(text) {
            return !!text && text.match(/^hsa/);
          });
        });
      })
      .reduce(function(accumulator, items) {
        items.forEach(function(item) {
          accumulator = accumulator.concat(item);
        });
        return accumulator;
      }, []);
    updateWidget(identifier, displayNames);
  });
