document.addEventListener('DOMContentLoaded', function() {
  const tocToggle = document.getElementById('tocToggle');
  const tocOverlay = document.getElementById('tocOverlay');
  const tocWrapper = document.querySelector('.tocify-wrapper') || document.querySelector('#TOC');

  if (tocToggle && tocWrapper) {
    // Toggle menu when clicking hamburger button
    tocToggle.addEventListener('click', function() {
      tocToggle.classList.toggle('active');
      tocWrapper.classList.toggle('active');
      tocOverlay.classList.toggle('active');
    });

    // Close menu when clicking overlay
    tocOverlay.addEventListener('click', function() {
      tocToggle.classList.remove('active');
      tocWrapper.classList.remove('active');
      tocOverlay.classList.remove('active');
    });

    // Close menu when clicking a TOC link
    const tocLinks = tocWrapper.querySelectorAll('a');
    tocLinks.forEach(link => {
      link.addEventListener('click', function() {
        tocToggle.classList.remove('active');
        tocWrapper.classList.remove('active');
        tocOverlay.classList.remove('active');
      });
    });

    // Close menu when pressing Escape key
    document.addEventListener('keydown', function(e) {
      if (e.key === 'Escape') {
        tocToggle.classList.remove('active');
        tocWrapper.classList.remove('active');
        tocOverlay.classList.remove('active');
      }
    });

    // Handle window resize - close menu if switching to large screen
    window.addEventListener('resize', function() {
      if (window.innerWidth > 1200) {
        tocToggle.classList.remove('active');
        tocWrapper.classList.remove('active');
        tocOverlay.classList.remove('active');
      }
    });
  }
});
